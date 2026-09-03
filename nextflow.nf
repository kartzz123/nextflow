nextflow.enable.dsl=2

// ============================================================
// PARAMETERS
// ============================================================

params.read1     = "/home/kartzz/ERR15695300_1.fastq"
params.read2     = "/home/kartzz/ERR15695300_2.fastq"
params.reference = "/home/kartzz/GENOME.fa"
params.outdir    = "/home/kartzz/results_out"

// ============================================================
// CLEAN / VALIDATE PAIRED-END FASTQ
// ============================================================

process CLEAN_FASTQ {

    tag "${sample_id}"

    publishDir "${params.outdir}/cleaned",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_1.clean.fastq.gz"),
          path("${sample_id}_2.clean.fastq.gz"),
          emit: reads

    script:
    """
    set -euo pipefail

    python3 <<'PY'
import gzip
import sys

r1_in  = "${r1}"
r2_in  = "${r2}"

r1_out = "${sample_id}_1.clean.fastq.gz"
r2_out = "${sample_id}_2.clean.fastq.gz"


def read_record(handle):
    h = handle.readline()

    if not h:
        return None

    seq  = handle.readline()
    plus = handle.readline()
    qual = handle.readline()

    # Detect truncated FASTQ records
    if not seq or not plus or not qual:
        return None

    return (
        h.rstrip("\\n\\r"),
        seq.rstrip("\\n\\r"),
        plus.rstrip("\\n\\r"),
        qual.rstrip("\\n\\r")
    )


def valid_record(record):
    if record is None:
        return False

    h, seq, plus, qual = record

    if not h.startswith("@"):
        return False

    if not plus.startswith("+"):
        return False

    if len(seq) == 0:
        return False

    if len(seq) != len(qual):
        return False

    return True


total = 0
kept = 0
removed = 0

with open(r1_in, "rt") as f1, \
     open(r2_in, "rt") as f2, \
     gzip.open(r1_out, "wt") as o1, \
     gzip.open(r2_out, "wt") as o2:

    while True:
        rec1 = read_record(f1)
        rec2 = read_record(f2)

        if rec1 is None and rec2 is None:
            break

        total += 1

        if not valid_record(rec1) or not valid_record(rec2):
            removed += 1
            continue

        o1.write(f"{rec1[0]}\\n{rec1[1]}\\n{rec1[2]}\\n{rec1[3]}\\n")
        o2.write(f"{rec2[0]}\\n{rec2[1]}\\n{rec2[2]}\\n{rec2[3]}\\n")
        kept += 1

print(f"Total pairs evaluated: {total}")
print(f"Pairs kept: {kept}")
print(f"Pairs discarded: {removed}")
PY
    """
}

// ============================================================
// INDEX REFERENCE GENOME
// ============================================================

process INDEX_REF {
    tag "${ref.baseName}"
    
    input:
    path ref

    output:
    path "${ref}*", emit: index

    script:
    """
    bwa index ${ref}
    # samtools command 1: Generate reference genome fasta index (.fai)
    samtools faidx ${ref}
    """
}

// ============================================================
// ALIGN READS AND CONVERT TO SORTED BAM
// ============================================================

process ALIGN_BWA {
    tag "${sample_id}"

    publishDir "${params.outdir}/alignment",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)
    path ref
    path ref_index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam_pair

    script:
    """
    set -euo pipefail
    # samtools command 2 & 3: Convert SAM stream to BAM view (-b) and pipe into sort
    bwa mem ${ref} ${r1} ${r2} | samtools view -bS - | samtools sort -o ${sample_id}.sorted.bam -
    
    # samtools command 4: Index the sorted BAM file (.bam.bai)
    samtools index ${sample_id}.sorted.bam
    """
}

// ============================================================
// SAMTOOLS QUALITY CONTROL AND ALIGNMENT STATISTICS
// ============================================================

process BAM_STATISTICS {
    tag "${sample_id}"

    publishDir "${params.outdir}/alignment_stats",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.flagstat.txt"),
          path("${sample_id}.stats.txt"),
          path("${sample_id}.idxstats.txt"),
          emit: qc_reports

    script:
    """
    set -euo pipefail
    
    # samtools command 5: Calculate quick alignment flag statistics
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    
    # samtools command 6: Run comprehensive metrics assessment (coverage, length, insert sizes)
    samtools stats ${bam} > ${sample_id}.stats.txt
    
    # samtools command 7: Report mapped/unmapped reads per reference chromosome/contig
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt
    """
}

// ============================================================
// VARIANT CALLING AND NORMALISATION
// ============================================================

process CALL_VARIANTS {
    tag "${sample_id}"

    publishDir "${params.outdir}/variants",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_index

    output:
    tuple val(sample_id), path("${sample_id}.norm.vcf.gz"), path("${sample_id}.norm.vcf.gz.tbi"), emit: vcf_pair

    script:
    """
    set -euo pipefail
    
    # Pileup and call variants
    bcftools mpileup -f ${ref} ${bam} | \
        bcftools call -mv -Ob -o ${sample_id}.raw.bcf

    # Left-align and normalize indels, then index the output file
    bcftools norm -f ${ref} -Oz -o ${sample_id}.norm.vcf.gz ${sample_id}.raw.bcf
    bcftools index ${sample_id}.norm.vcf.gz
    """
}

// ============================================================
// GENERATE CONSENSUS FASTA
// ============================================================

process GENERATE_CONSENSUS {
    tag "${sample_id}"

    publishDir "${params.outdir}/consensus",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path ref
    path ref_index

    output:
    path "${sample_id}_consensus.fa", emit: consensus

    script:
    """
    set -euo pipefail
    
    # Apply normalized variants to the reference genome to create consensus sequence
    bcftools consensus -f ${ref} ${vcf} > ${sample_id}_consensus.fa
    """
}

// ============================================================
// MAIN WORKFLOW ENTRY POINT
// ============================================================

workflow {
    def sample_id = file(params.read1).baseName.replaceAll(/_1$/, "")
    
    reads_ch = Channel.of( tuple(sample_id, file(params.read1), file(params.read2)) )
    ref_ch   = Channel.fromPath(params.reference, checkIfExists: true)

    // Step 1: Clean Reads
    CLEAN_FASTQ(reads_ch)

    // Step 2: Index reference genome (generates BWA indexes & .fai index)
    INDEX_REF(ref_ch)

    // Step 3: Align, convert to BAM, sort, and index
    ALIGN_BWA(CLEAN_FASTQ.out.reads, ref_ch, INDEX_REF.out.index.first())

    // Step 4: Generate Samtools alignment metrics and statistics reports
    BAM_STATISTICS(ALIGN_BWA.out.bam_pair)

    // Step 5: Run bcftools mpileup, call, and norm
    CALL_VARIANTS(ALIGN_BWA.out.bam_pair, ref_ch, INDEX_REF.out.index.first())

    // Step 6: Run bcftools consensus to generate consensus.fa
    GENERATE_CONSENSUS(CALL_VARIANTS.out.vcf_pair, ref_ch, INDEX_REF.out.index.first())
}
