#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 30:00
#SBATCH -J quast_scaffold
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/11_scaffolding/04_evaluation/%x.%j.out

module load bioinfo-tools
module load quast/5.0.2

SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="$SRCDIR/analyses/11_scaffolding/04_evaluation/"
SCAFFOLDS="$SRC_DIR/analyses/11_scaffolding/02_scaffolding/yahs.out_scaffolds_final.fa"

cd $JOB_DIR

python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py $SCAFFOLDS \
  --eukaryote \
  --threads 4 \
  --large \
  --labels Scaffold \
  -o $JOB_DIR

