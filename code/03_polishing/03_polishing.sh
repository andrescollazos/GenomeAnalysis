#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-00:00:00
#SBATCH -J anco_polishing
#SBATCH --mail-type=ALL
#SBATCH --mail-user andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/03_polishing/%x.%j.out

module load bioinfo-tools
module load samtools/1.20
module load bwa/0.7.18
module load Pilon/1.24

export SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
export JOB_DIR="$SRCDIR/analyses/03_polishing/"

export REF_GENOME="assembly.fasta"
export FORWARD_READ="$SRCDIR/data/raw/dna_reads/chr3_illumina_R1.fastq.gz"
export REVERSE_READ="$SRCDIR/data/raw/dna_reads/chr3_illumina_R2.fastq.gz"


# Creates an allignment for the illumina reads using the assembled genome as a reference
cd $JOB_DIR/01_mapping
bwa index $REF_GENOME
bwa mem -t 4 $REF_GENOME $FORWARD_READ $REVERSE_READ | samtools view -Sb - | samtools sort -o$
samtools index illumina_sorted.bam

# POLISHING
cd ..
mkdir -p 02_polishing
cd 02_polishing

java -Xmx16G -jar $PILON_HOME/pilon.jar \
  --genome $REF_GENOME \
  --frags ../01_mapping/illumina_sorted.bam \
  --output polished_assembly \
  --outdir pilon_output \
  --threads 4
