#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 3:00:00
#SBATCH -J masking_with_dfam_lib
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/05_masking/03_masking_with_species/%x.%j.out

module load bioinfo-tools
module load RepeatMasker/4.1.5

export SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/05_masking/03_masking_with_species"
export ASSEMBLY="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

mkdir -p "$JOB_DIR"
cd "$JOB_DIR"

RepeatMasker \
  -species "Bryophyta" \
  -pa 8 \
  -xsmall \
  -no_is \
  -gff \
  -html \
  -dir $JOB_DIR \
  "$ASSEMBLY"

# -species "Bryophyta"    From Dfam
# -pa 8                 8 threads 
# -xsmall               Mask repeats in lowercase 
# -no_is                Skip bacterial insertion elements check
# -gff                  Output GFF3 format annotation file

