#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J eggNOGMapper
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/08_functional_annotation/02_eggNOG_mapper/%x.%j.out

# Load Modules
module load bioinfo-tools eggNOG-mapper/2.1.9

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/08_functional_annotation/02_eggNOG_mapper/"

PROTEINS_FILE="$SRC_DIR/analyses/07_structural_annotation/braker_standardized.aa"
GENOME_FILE="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

emapper.py -i "$PROTEINS_FILE" \
  --itype proteins \
  -m diamond \
  --cpu 8 \
  --go_evidence experimental \
  --output n_japonicum \
  --output_dir "$JOB_DIR" \
  --decorate_gff "$SRC_DIR/analyses/07_structural_annotation/braker_standardized.gff3" \
  --decorate_gff_ID_field ID \
  --excel