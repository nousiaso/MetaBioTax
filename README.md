# MetaBioTax
MetaBioTax is a WGS (whole genome shotgun) metagenomics pipeline by Orestis Nousias (with eukaryotes in mind), designed to align Illumina and Nanopore reads, classify sequences taxonomically, and analyze biodiversity. It efficiently processes genomic data, generating taxonomic summaries and species counts for key taxa, optimized for large-scale metagenomic studies.

Overview

MetaBioTax is a comprehensive metagenomics pipeline designed to perform sequence alignment, taxonomic classification, and biodiversity analysis of both Illumina and Nanopore sequencing data. This pipeline efficiently processes genomic data, annotates taxonomic diversity, and generates detailed species-level summaries, offering an all-in-one solution for metagenomic and taxonomic biodiversity research.

Author: Orestis Nousias

Features

Supports Multiple Sequencing Types: Aligns both paired-end Illumina reads and single-end Nanopore reads against the NR protein database.
Taxonomic Classification: Annotates sequences with taxonomic data using the NCBI taxonomy database, providing a hierarchical breakdown of biodiversity.
Species-Specific Analysis: Generates counts for species of interest (e.g., mammals and metazoans) from aligned sequences.
Customizable Parameters: Allows fine-tuning of alignment settings for enhanced sensitivity and specificity.
Parallel Processing: Optimized to run on multi-core systems for fast and efficient data processing.
Pipeline Workflow

File Downloads: The pipeline automatically downloads the necessary NCBI databases:
nr.gz: NR protein database.
nodes.dmp and names.dmp: NCBI taxonomy nodes and names files.
prot.accession2taxid: Protein accession to taxonomic ID mapping file.
Data Alignment:
Paired-End Illumina Reads: Reads are aligned against the NR database using DIAMOND's blastx mode, with taxonomic annotations generated.
Nanopore Single-End Reads: Nanopore sequences are aligned with additional parameters (-F 15, --range-culling, --top 10) to capture top hits and manage sequence range culling.
Taxonomic and Species Analysis:
Taxonomic Frequency: A Python subscript processes the BLAST results to generate a summary of the taxonomic counts for specified taxa.
Mammal and Metazoan Species Counting: Separate scripts generate counts for mammal and metazoan species, which are output to respective files.
Input Requirements

Paired-End Illumina FASTQ files: For example, R1.fastq.gz and R2.fastq.gz.
Single-End Nanopore FASTQ file: Example, nanopore.fastq.gz.
Species Lists: Plain text files containing species of interest for biodiversity analysis (one species per line).
Output

Aligned Data: DIAMOND alignment files in DAA format.
Taxonomic Summaries: A summary of taxonomic data for specific taxa in summary.txt.
Species Counts: Mammal and metazoan species counts in mammal_counts.txt and metazoa_count.txt, respectively.
System Requirements

DIAMOND v2.1.8 (or later).
Python 3 with required libraries (os, re, collections).
Multi-core system (recommended: 20 CPU cores and 40GB RAM).
Access to high-performance computing for large datasets.

Clone the repository and ensure all dependencies are installed.
Modify the SLURM directives as needed for your system.
Provide the necessary input files (FASTQ reads and species lists).

Check the Supplementary table for existing fastq ssequences. These can be used for testing purposes of the pipeline. Running thie pipeline on a "normal" pc, is not advised. It will take multiple days to finish on 10 cores.
