#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 30:00
#SBATCH -J anco_pre_quality_control
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/01_preprocessing/01_before_trimming/%x.%j.out

module load bioinfo-tools
module load FastQC/0.11.2

export SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export OUTPUT_DIR="$SRCDIR/analyses/01_preprocessing/01_before_trimming"
export FORWARD="$SRCDIR/data/raw/dna_reads/chr3_illumina_R1.fastq.gz"
export REVERSE="$SRCDIR/data/raw/dna_reads/chr3_illumina_R2.fastq.gz"

fastqc $FORWARD $REVERSE --outdir=$OUTPUT_DIR -t 2

