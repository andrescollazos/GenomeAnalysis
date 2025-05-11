#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -J rna_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/06_rna_mapping/01_mapping/%x.%j.out

# Load Modules
module load bioinfo-tools
module load HISAT2/2.2.1

export SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/06_rna_mapping/01_mapping"
export ASSEMBLY="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"
export RNA_READS="$SRC_DIR/data/rna_reads/"

mkdir -p $JOB_DIR
cd "$JOB_DIR"

# Build the Index
hisat2-build "$ASSEMBLY" chr3_index


for condition in "Control" "Heat_treated_42_12h"
do
 for replicate in 1 2 3
  do
   echo Processing "$condition"_"$replicate"_f1.fq.gz and "$condition"_"$replicate"_r2.fq.gz
   hisat2 --threads 8 -x chr3_index \
	-q \
	-1 "$RNA_READS/${condition}_${replicate}_f1.fq.gz" \
	-2 "$RNA_READS/${condition}_${replicate}_r2.fq.gz" \
	-S "${condition}_${replicate}.sam"
	2> "${condition}_${replicate}.log"
  done
done

