# Streamlined chr20-Based Variant-Calling Pipeline

A streamlined, chr20-focused variant-calling pipeline built using real Illumina sequencing data (SRR1910621).
Designed to demonstrate all major steps of an NGS variant-calling workflow while remaining computationally efficient and suitable for systems with modest resources (e.g., 8 GB RAM).

# Project Overview

This workflow performs:

FASTQ acquisition

Quality control and trimming

Alignment to a chr20-only reference

Duplicate marking

Variant calling (bcftools)

Variant filtering

Variant statistics

Variant density analysis

Restricting analysis to chromosome 20 enables rapid experimentation while preserving biological interpretability.

# Directory Structure
ngs_variant_project/
├── raw_data/                 # FASTQ files (excluded from GitHub)
├── reference/
│   └── chr20.sizes           # Reference size file (included)
├── qc/
│   └── fastp_report.html     # QC report
├── alignment/                # BAM + index files (excluded)
├── variants/                 # VCF files (excluded)
├── analysis/
│   ├── variant_stats.txt
│   ├── chr20_100kb_windows.bed
│   └── variant_density_100kb.bed
└── .gitignore



Large FASTQs, BAMs, VCFs, and reference FASTA/index files are intentionally excluded with .gitignore.

# Dataset

Accession: SRR1910621
Sequencing: Illumina paired-end
Source: NCBI SRA
Size: ~600 MB compressed

# Tools Used

fastp – QC and trimming

bwa mem – alignment

samtools – sorting & indexing

Picard – duplicate marking

bcftools – variant calling & filtering

bedtools – windowing & density analysis

# Pipeline Summary
## 1. Quality Control
fastp \
  -i raw_data/sample_R1.fastq.gz \
  -I raw_data/sample_R2.fastq.gz \
  -o qc/clean_R1.fastq.gz \
  -O qc/clean_R2.fastq.gz \
  -h qc/fastp_report.html

## 2. Alignment to chr20
bwa mem -t 4 \
  -R '@RG\tID:1\tSM:sample\tLB:lib1\tPL:ILLUMINA' \
  reference/chr20.fa \
  qc/clean_R1.fastq.gz qc/clean_R2.fastq.gz \
| samtools sort -o alignment/sample.chr20.sorted.bam

samtools index alignment/sample.chr20.sorted.bam

## 3. Mark Duplicates
picard MarkDuplicates \
  I=alignment/sample.chr20.sorted.bam \
  O=alignment/sample.chr20.dedup.bam \
  M=alignment/dup_metrics.txt

samtools index alignment/sample.chr20.dedup.bam

## 4. Variant Calling
bcftools mpileup \
  -f reference/chr20.fa \
  -a FORMAT/AD \
  alignment/sample.chr20.dedup.bam \
| bcftools call -mv -Ov \
  -o variants/sample.chr20.raw.vcf

## 5. Variant Filtering
bcftools filter \
  -s LOWQUAL \
  -e 'DP<10 || QUAL<20' \
  variants/sample.chr20.raw.vcf \
  -Ov -o variants/sample.chr20.filtered.vcf

## 6. Variant Statistics
bcftools stats variants/sample.chr20.filtered.vcf \
  > analysis/variant_stats.txt

## 7. Variant Density
bedtools makewindows -g reference/chr20.sizes -w 100000 \
  > analysis/chr20_100kb_windows.bed

bedtools intersect \
  -a analysis/chr20_100kb_windows.bed \
  -b variants/sample.chr20.filtered.vcf \
  -c > analysis/variant_density_100kb.bed

# Included Outputs

fastp_report.html – QC summary

variant_stats.txt – SNP/INDEL counts & Ti/Tv ratio

variant_density_100kb.bed – variant distribution

chr20_100kb_windows.bed – coordinate windows

# Reproducibility

To reproduce this workflow:

Clone this repository

Download chr20.fa (hg38)

Download SRR1910621 FASTQs

Run each pipeline command in order

# Future Extensions

R-based visualization

Ti/Tv and base-change summary plots

Depth & coverage plots

Gene-level annotation with SnpEff or VEP

Multi-sample comparisons
