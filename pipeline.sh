#!/bin/bash

#SBATCH --job-name=aWL1wJUL24               # Job name for tracking purposes
#SBATCH --output=aWL1wJUL24.txt             # Output file to capture stdout
#SBATCH --time=10-00:00:00                  # Max time for job execution (10 days)
#SBATCH --mail-type=ALL                     # Notifications for all job events
#SBATCH --ntasks=1                          # Number of tasks to run
#SBATCH --cpus-per-task=20                  # Number of CPU cores allocated
#SBATCH --mem=40G                           # Memory allocation (40GB)

# Load necessary modules
module load diamond/2.1.8                   # Load DIAMOND tool for protein alignment

# Step 1: Download necessary files for NR database and taxonomy mapping

# Download nr.gz (NR database for protein sequences)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -O nr.gz
# Download nodes.dmp (NCBI taxonomy nodes file)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O taxdump.tar.gz
tar -xvzf taxdump.tar.gz nodes.dmp
# Download names.dmp (NCBI taxonomy names file)
tar -xvzf taxdump.tar.gz names.dmp
# Download prot.accession2taxid (Protein accession to taxonomic ID mapping)
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz -O prot.accession2taxid.gz
gunzip prot.accession2taxid.gz

# Set directory and file variables
READ1="aWL1wJul24_S11_L008_R1_001.fastq.gz" # First read file (FASTQ format, paired-end)
READ2="aWL1wJul24_S11_L008_R2_001.fastq.gz" # Second read file (FASTQ format, paired-end)
NANOPORE_READ="aWL1wJul24_nano.fastq.gz"    # Nanopore single-read file (FASTQ format)
DB_NR="nr.gz"                               # NR database file
PROT_ACC2TAXID="prot.accession2taxid"       # Protein accession to taxonomic ID file
NODES="nodes.dmp"                           # Taxonomic nodes file
NAMES="names.dmp"                           # Taxonomic names file
THREADS=20                                  # Number of threads to use for parallelization

# Step 2: Process data using minimap2 (optional, commented out for now)

## 2.1. Indexing the NR database using minimap2 (commented out as DIAMOND is preferred here)
#minimap2 -d nr.mmi $DB_NR

## 2.2. Align reads using minimap2 (commented out)
#minimap2 -ax sr -t $THREADS nr.mmi $READ1 $READ2 > aligned_to_nr.sam

# Step 3: Process data using DIAMOND for NR database

## 3.1. Build DIAMOND database with taxonomic information for NR
diamond makedb --in $DB_NR -d nr --taxonmap $PROT_ACC2TAXID --taxonnodes $NODES --taxonnames $NAMES

## 3.2. Align paired-end reads (Illumina) against NR database with taxonomic information
diamond blastx -d nr -q "$READ1" "$READ2" -o aligned_to_nr_awL1wJul24.daa -p $THREADS --outfmt 6 \
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle --max-target-seqs 1 --more-sensitive

## 3.3. Align nanopore single read sequences against NR database with additional parameters
diamond blastx -d nr -q "$NANOPORE_READ" -o aligned_to_nr_nanopore.daa -p $THREADS --outfmt 6 \
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle --max-target-seqs 10 --more-sensitive \
-F 15 --range-culling --top 10

# Completion message after successful alignment
echo "DIAMOND alignment against NR database is complete for both Illumina and Nanopore reads."

# Subscript 1: Python script to process taxonomic BLAST output for specific taxa
python3 - << 'END_SCRIPT'

import os

def get_tax_id_from_names_dmp(taxon_name, names_dmp_path):
    print(f"Searching for taxon: {taxon_name}")
    with open(names_dmp_path, 'r') as f:
        for line in f:
            fields = line.split('|')
            if fields[1].strip() == taxon_name:
                return fields[0].strip()
    return None

def read_blast_output(blast_output_file):
    print("Reading BLAST output")
    tax_id_freq = {}
    with open(blast_output_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            staxid = fields[12]
            tax_id_freq[staxid] = tax_id_freq.get(staxid, 0) + 1
    return tax_id_freq

def build_parent_dict(nodes_dmp_path):
    print("Building parent dictionary")
    parent_dict = {}
    with open(nodes_dmp_path, 'r') as f:
        for line in f:
            fields = line.split('|')
            tax_id = fields[0].strip()
            parent_id = fields[1].strip()
            parent_dict[tax_id] = parent_id
    return parent_dict

# Initialize
taxa_list = ["Arthropoda", "Mollusca", "Annelida", "Nematoda", "Platyhelminthes", "Rotifera", 
             "Gastrotricha", "Bryozoa", "Brachiopoda", "Phoronida", "Sipuncula", "Priapulida", 
             "Tardigrada", "Chordata", "Echinodermata", "Hemichordata", "Xenoturbellida", 
             "Acoelomorpha"]
names_dmp_path = os.path.join(os.getcwd(), "names.dmp")
nodes_dmp_path = os.path.join(os.getcwd(), "nodes.dmp")
blast_output_file = os.path.join(os.getcwd(), 'aligned_to_nr_awL1wJul24.daa')

# Main Execution
tax_ids_of_interest = {taxon: get_tax_id_from_names_dmp(taxon, names_dmp_path) for taxon in taxa_list}
tax_id_freq = read_blast_output(blast_output_file)
parent_dict = build_parent_dict(nodes_dmp_path)
summary_dict = {taxon: 0 for taxon in taxa_list}

for staxid, freq in tax_id_freq.items():
    print(f"Processing staxid: {staxid}, freq: {freq}")
    current_tax_id = staxid
    max_iter = 1000000  # Adjust based on a reasonable expectation
    iter_count = 0

    while current_tax_id in parent_dict:
        iter_count += 1
        if iter_count > max_iter:
            print(f"Exiting to prevent infinite loop for staxid: {staxid}")
            break
        if current_tax_id in tax_ids_of_interest.values():
            taxon = [k for k, v in tax_ids_of_interest.items() if v == current_tax_id][0]
            summary_dict[taxon] += freq
            print(f"Incremented count for {taxon} by {freq}")
            break
        current_tax_id = parent_dict[current_tax_id]

with open('summary_awL1wJul24.txt', 'w') as f:
    for taxon, freq in summary_dict.items():
        f.write(f"{taxon}\t{freq}\n")

print("Taxa summary has been written to 'summary_awL1wJul24.txt'")
END_SCRIPT

# Subscript 2: Python script to count mammal species in BLAST output
python3 - << 'END_SCRIPT'

import os
import re
from collections import Counter

def read_mammal_species(file_path):
    with open(file_path, 'r') as file:
        return set(line.strip() for line in file)

def read_blast_output(blast_output_file, species_set):
    print("Reading BLAST output")
    species_freq = Counter()
    pattern = re.compile(r'\[(.*?)\]')  # Regex pattern for extracting species name
    with open(blast_output_file, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                species_name = match.group(1).strip()
                if species_name in species_set:
                    species_freq[species_name] += 1
    return species_freq

# Initialize
mammal_species_list_path = os.path.join(os.getcwd(), "mammal_species_list.txt")
mammal_species_set = read_mammal_species(mammal_species_list_path)
blast_output_file = os.path.join(os.getcwd(), 'aligned_to_nr_awL1wJul24.daa')

# Main Execution
species_freq = read_blast_output(blast_output_file, mammal_species_set)

# Sort species by frequency, highest to lowest, and filter out zero counts
sorted_species = sorted(species_freq.items(), key=lambda x: x[1], reverse=True)
non_zero_species = [(species, count) for species, count in sorted_species if count > 0]

# Write results to file
with open('awL1wJul24_mammal_counts.txt', 'w') as f:
    for species, freq in non_zero_species:
        f.write(f"{species}\t{freq}\n")

print("Mammal species count has been written to 'awL1wJul24_mammal_counts.txt'")
END_SCRIPT

# Subscript 3: Python script to count metazoan species in BLAST output
python3 - << 'END_SCRIPT'

import os
import re
from collections import Counter

def read_metazoa_species(file_path):
    with open(file_path, 'r') as file:
        return set(line.strip() for line in file)

def read_blast_output(blast_output_file, species_set):
    print("Reading BLAST output")
    species_freq = Counter()
    pattern = re.compile(r'\[(.*?)\]')  # Regex pattern for extracting species name
    with open(blast_output_file, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                species_name = match.group(1).strip()
                if species_name in species_set:
                    species_freq[species_name] += 1
    return species_freq

# Initialize
metazoa_species_list_path = os.path.join(os.getcwd(), "metazoa_species_list.txt")
metazoa_species_set = read_metazoa_species(metazoa_species_list_path)
blast_output_file = os.path.join(os.getcwd(), 'aligned_to_nr_awL1wJul24.daa')

# Main Execution
species_freq = read_blast_output(blast_output_file, metazoa_species_set)

# Sort species by frequency, highest to lowest, and filter out zero counts
sorted_species = sorted(species_freq.items(), key=lambda x: x[1], reverse=True)
non_zero_species = [(species, count) for species, count in sorted_species if count > 0]

# Write results to file
with open('awL1wJul24_metazoa_count.txt', 'w') as f:
    for species, freq in non_zero_species:
        f.write(f"{species}\t{freq}\n")

print("Metazoan species count has been written to 'awL1wJul24_metazoa_count.txt'")
END_SCRIPT
