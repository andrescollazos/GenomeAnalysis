#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -J custom_masking_library
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/05_masking/%x.%j.out

module load bioinfo-tools
module load RepeatMasker/4.1.5

export SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="$SRC_DIR/analyses/05_masking"
export ASSEMBLY="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

mkdir -p $JOB_DIR
cd $JOB_DIR

RepeatMasker \
  -species "Bryophyta" \  # From Dfam
  -pa 8 \                 # 8 threads 
  -xsmall \               # Mask repeats in lowercase 
  -no_is \                # Skip bacterial insertion elements check
  -gff \                  # Output GFF3 format annotation file
  -dir $JOB_DIR \          
  $ASSEMBLY

