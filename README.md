Paired-End Read Alignment, Variant Calling & Consensus Genome Pipeline

A Nextflow DSL2 bioinformatics pipeline for processing paired-end FASTQ sequencing data, aligning reads to a reference genome, performing alignment quality control, calling and normalizing variants, and generating a sample-specific consensus FASTA sequence.

The workflow integrates BWA, Samtools, and BCFtools into a reproducible, modular analysis pipeline.

Workflow Overview

The pipeline performs the following steps:

Paired-End FASTQ
       │
       ▼
┌─────────────────────┐
│   FASTQ Cleaning    │
│  Validate read pairs│
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐       ┌─────────────────────┐
│    BWA Alignment    │◄──────│ Reference Indexing  │
│   bwa mem + sort    │       │ BWA index + .fai    │
└──────────┬──────────┘       └─────────────────────┘
           │
           ▼
┌─────────────────────┐
│ Alignment QC/Stats  │
│ flagstat/stats/     │
│ idxstats             │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│   Variant Calling   │
│ mpileup + call +    │
│ normalization        │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Consensus Generation │
│ bcftools consensus  │
└──────────┬──────────┘
           │
           ▼
    Consensus FASTA

Pipeline Steps
1. FASTQ Cleaning and Validation

The CLEAN_FASTQ process reads paired-end FASTQ files and validates individual FASTQ records.

It checks:

FASTQ header starts with @
+ separator line is present
Sequence is non-empty
Sequence and quality strings have equal lengths
Both mates contain valid records

Invalid read pairs are discarded, while valid pairs are written as compressed .fastq.gz files.

The process also reports:

Total pairs evaluated
Pairs kept
Pairs discarded

2. Reference Genome Indexing

The INDEX_REF process prepares the reference genome for downstream analysis using:

bwa index — creates BWA alignment indexes
samtools faidx — creates the FASTA index (.fai)
3. Read Alignment

The ALIGN_BWA process aligns the cleaned paired-end reads to the reference genome using bwa mem.

The SAM output is streamed directly into Samtools to:

Convert SAM to BAM
Sort the BAM file
Create a BAM index

The resulting files are:

sample.sorted.bam
sample.sorted.bam.bai

4. Alignment Quality Control

The BAM_STATISTICS process generates three Samtools reports:

sample.flagstat.txt
sample.stats.txt
sample.idxstats.txt


These provide information about:

Total and mapped reads
Alignment statistics
Coverage-related metrics
Insert-size statistics
Reference/contig-level mapping statistics
Mapped and unmapped reads
5. Variant Calling and Normalization

The CALL_VARIANTS process uses BCFtools to identify variants from the aligned BAM file.

The workflow performs:

bcftools mpileup
        ↓
bcftools call
        ↓
bcftools norm


Variants are initially written to BCF and subsequently normalized and compressed into:

sample.norm.vcf.gz
sample.norm.vcf.gz.tbi


Variant normalization helps standardize representation of variants, particularly indels, relative to the reference genome.

6. Consensus Genome Generation

The GENERATE_CONSENSUS process applies the called variants to the reference genome using:

bcftools consensus


The final output is:

sample_consensus.fa


This represents the sample's consensus sequence based on the reference genome and the variants identified from the sequencing data.

Output Directory Structure

Results are organized automatically using Nextflow's publishDir:

results_out/
├── cleaned/
│   ├── sample_1.clean.fastq.gz
│   └── sample_2.clean.fastq.gz
│
├── alignment/
│   ├── sample.sorted.bam
│   └── sample.sorted.bam.bai
│
├── alignment_stats/
│   ├── sample.flagstat.txt
│   ├── sample.stats.txt
│   └── sample.idxstats.txt
│
├── variants/
│   ├── sample.norm.vcf.gz
│   └── sample.norm.vcf.gz.tbi
│
└── consensus/
    └── sample_consensus.fa

Software Requirements

The pipeline requires:

Nextflow
Python 3
BWA
Samtools
BCFtools

Make sure all tools are available in the execution environment and accessible through $PATH.

Input Parameters

The pipeline currently defines the following parameters:

params.read1     = "/home/kartzz/ERR15695300_1.fastq"
params.read2     = "/home/kartzz/ERR15695300_2.fastq"
params.reference = "/home/kartzz/GENOME.fa"
params.outdir    = "/home/kartzz/results_out"


These can be modified directly in the workflow or overridden from the command line.

For example:

nextflow run main.nf \
    --read1 /path/to/sample_1.fastq \
    --read2 /path/to/sample_2.fastq \
    --reference /path/to/reference.fa \
    --outdir /path/to/results

Running the Pipeline

Clone the repository and execute:

nextflow run main.nf


To resume an interrupted or failed workflow:

nextflow run main.nf -resume

Key Features
Nextflow DSL2 modular workflow architecture
Paired-end FASTQ validation and cleaning
Automated reference genome indexing
BWA-based read alignment
Sorted and indexed BAM generation
Samtools alignment quality-control reports
BCFtools variant calling
Variant normalization
Sample-specific consensus FASTA generation
Organized output directories
Support for Nextflow's caching and -resume functionality
Workflow Components
Process	Purpose	Main Tools
CLEAN_FASTQ	Validate and compress paired FASTQ reads	Python
INDEX_REF	Index reference genome	BWA, Samtools
ALIGN_BWA	Align reads and generate sorted BAM	BWA, Samtools
BAM_STATISTICS	Generate alignment QC statistics	Samtools
CALL_VARIANTS	Call and normalize variants	BCFtools
GENERATE_CONSENSUS	Generate consensus genome	BCFtools
Intended Use

This pipeline is designed as a compact workflow for short-read resequencing analysis, particularly applications where the goal is to progress from raw paired-end FASTQ reads to:

FASTQ → Clean Reads → BAM → Variants → Consensus FASTA


It can serve as a starting point for human genome analysis, targeted resequencing, or other reference-based variant analysis workflows.

Important Notes

This pipeline is intentionally lightweight and does not currently include adapter trimming, quality-score-based read trimming, duplicate marking, base-quality recalibration, variant filtering, or depth/coverage-based variant filtering.

For production analyses, additional validation and filtering steps may be appropriate depending on the organism, sequencing platform, experimental design, and downstream objectives.

License

Add an appropriate license for your project, such as MIT, GPL-3.0, or Apache-2.0.
