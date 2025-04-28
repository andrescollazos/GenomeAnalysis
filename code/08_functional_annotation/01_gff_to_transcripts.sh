#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J gff3_to_protein
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/08_functional_annotation/01_get_transcripts/%x.%j.out

# Load Modules
module load bioinfo-tools gffread/0.12.7

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="$SRC_DIR/analyses/08_functional_annotation/01_get_transcripts/"

# Set input files
GFF3_FILE="$SRC_DIR/analyses/07_structural_annotation/braker.gff3"
GENOME_FILE="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

# Extract protein sequences
gffread -y "${JOB_DIR}/braker.aa" -g "$GENOME_FILE" "$GFF3_FILE" 
