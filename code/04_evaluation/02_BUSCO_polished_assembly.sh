#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -J busco_polished
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/04_evaluation/02_BUSCO_polished_assembly/%x.%j.out

module load bioinfo-tools
module load BUSCO/5.7.1

# Paths
export SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="$SRCDIR/analyses/04_evaluation/02_BUSCO_polished_assembly"
export ASSEMBLY="$JOB_DIR/polished_assembly.fasta"
cd $JOB_DIR

# Set up AUGUSTUS config locally
source $AUGUSTUS_CONFIG_COPY

# Run BUSCO
busco -i $ASSEMBLY \
      -l $BUSCO_LINEAGE_SETS/embryophyta_odb10/ \
      -o busco_polished \
      -m genome \
      -c 4 \
      -f


# Embryophyta: Land Plants
