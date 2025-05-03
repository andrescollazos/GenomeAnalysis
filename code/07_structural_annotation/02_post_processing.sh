#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J post_struc_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/07_structural_annotation/%x.%j.out

# Load Modules
module load bioinfo-tools AGAT/1.3.2

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/07_structural_annotation/braker"

GFF3_FILE="$JOB_DIR/braker.gff3"
GFF_CLEANED="$JOB_DIR/braker_cleaned.gff3"
GFF_STD_NAME="braker_standardized.gff3"
GFF_STANDARD_PATH="$JOB_DIR/$GFF_STD_NAME"

GENOME_FILE="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

# convert GFF into GFF (bioperl) format
agat_convert_sp_gxf2gxf.pl -gff $GFF3_FILE -o $GFF_CLEANED

# manage IDs
agat_sp_manage_IDs.pl -gff $GFF_CLEANED -o $GFF_STANDARD_PATH

# Get feature statistics
agat_sp_statistics.pl --gff $GFF_STANDARD_PATH -o "$JOB_DIR/braker_standardized.statistics.txt"

# Extract protein sequences
agat_sp_extract_sequences.pl --gff $GFF_STANDARD_PATH --fasta $GENOME_FILE -p -o "$JOB_DIR/braker_standardized.aa"

# Create symlinks
OUT_DIR="$SRC_DIR/analyses/07_structural_annotation"
ln -sf $JOB_DIR/*_standardized.* $OUT_DIR
