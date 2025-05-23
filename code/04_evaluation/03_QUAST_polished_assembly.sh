#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 30:00
#SBATCH -J quast_polished
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/04_evaluation/03_QUAST_polished_assembly/%x.%j.out

module load bioinfo-tools
module load quast/5.0.2

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="$SRC_DIR/analyses/04_evaluation/03_QUAST_polished_assembly"
ASSEMBLY="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"

cd $JOB_DIR

python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py $ASSEMBLY \
  --eukaryote \
  --threads 4 \
  --large \
  --labels Polished \
  -o $JOB_DIR

