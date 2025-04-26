#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -J custom_masking_library
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/05_masking/01_customer_library/%x.%j.out

# Load Modules
module load bioinfo-tools
module load RepeatModeler/2.0.4

export SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="$SRC_DIR/analyses/05_masking/01_customer_library"
export ASSEMBLY="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

mkdir -p "$JOB_DIR"
cd "$JOB_DIR"

# Step 1: Build RepeatModeler database
BuildDatabase -name njaponicum_db -engine ncbi $ASSEMBLY

# Step 2: Run RepeatModeler with LTRStruct option
RepeatModeler -database njaponicum_db -threads 8 -LTRStruct
