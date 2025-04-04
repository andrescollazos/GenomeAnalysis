#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J anco_flye_assembly​
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se

module load bioinfo-tools Flye/2.9.5

export SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export NANO_FILE="chr3_clean_nanopore.fq.gz"
export INPUT_FILE="$SRCDIR/data/raw/dna_reads/$NANO_FILE"
export OUTPUT_DIR="$SRCDIR/analyses/02_genome_assembly/flye_out"

flye --nano-raw $INPUT_FILE --out-dir $OUTPUT_DIR --threads 8 --resume

