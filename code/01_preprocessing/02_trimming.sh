#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 45:00
#SBATCH -J trimming
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/01_preprocessing/02_trimming/%x.%j.out

module load bioinfo-tools
module load trimmomatic/0.39

export SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export OUTPUT_DIR="$SRCDIR/analyses/01_preprocessing/02_trimming"

export IN_FORWARD="$SRCDIR/data/raw/dna_reads/chr3_illumina_R1.fastq.gz"
export IN_REVERSE="$SRCDIR/data/raw/dna_reads/chr3_illumina_R2.fastq.gz"
export OUT_FORWARD_PAIRED="$OUTPUT_DIR/chr3_illumina_R1_paired.fastq.gz"
export OUT_FORWARD_UNPAIRED="$OUTPUT_DIR/chr3_illumina_R1_unpaired.fastq.gz"
export OUT_REVERSE_PAIRED="$OUTPUT_DIR/chr3_illumina_R2_paired.fastq.gz"
export OUT_REVERSE_UNPAIRED="$OUTPUT_DIR/chr3_illumina_R2_unpaired.fastq.gz"
export OUT_LOG="$OUTPUT_DIR/PE_trimmomatic.log"

trimmomatic PE -threads 2 -phred33 -trimlog $OUT_LOG \
    $IN_FORWARD $IN_REVERSE \
    $OUT_FORWARD_PAIRED $OUT_FORWARD_UNPAIRED $OUT_REVERSE_PAIRED $OUT_REVERSE_UNPAIRED \
    ILLUMINACLIP:$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40


