#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J assembly​
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/02_genome_assembly/flye_out/%x.%j.out

module load bioinfo-tools
module load Flye/2.95.

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
NANOPORE="$SRC_DIR/data/raw/dna_reads/chr3_clean_nanopore.fq.gz"
OUTPUT_DIR="$SRC_DIR/analyses/02_genome_assembly/flye_out"

flye --nano-raw $NANOPORE --out-dir $OUTPUT_DIR --threads 8 --resume