#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 30:00
#SBATCH -J quast_unpolished
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/04_evaluation/01_QUAST_unpolished_assembly/%x.%j.out

module load bioinfo-tools
module load quast/5.0.2

export SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="$SRCDIR/analyses/04_evaluation/01_QUAST_unpolished_assembly"
export RAW_ASSEMBLY="./assembly.fasta"

cd $JOB_DIR

python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py $RAW_ASSEMBLY \
  --eukaryote \
  --threads 2 \
  --large \
  --labels Unpolished \
  -o $JOB_DIR

