#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J hic_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/11_scaffolding/01_hic_mapping/%x.%j.out

module load bioinfo-tools
module load samtools/1.20
module load bwa/0.7.18

THREADs=4

SRCDIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/11_scaffolding/01_hic_mapping"

REF_GENOME="$JOB_DIR/polished_assembly.fasta"
FORWARD_READ="$SRCDIR/data/raw/hi_reads/chr3_hiC_R1.fastq.gz"
REVERSE_READ="$SRCDIR/data/raw/hi_reads/chr3_hiC_R2.fastq.gz"

cd $JOB_DIR

bwa index $REF_GENOME
bwa mem -t $THREADs $REF_GENOME $FORWARD_READ $REVERSE_READ | \
samtools view -Sb - > aln.bam

samtools collate -@ $THREADs -O -u aln.bam | \
samtools fixmate -@ $THREADs -m -u - - | \
samtools sort -@ $THREADs -u - | \
samtools markdup -@ $THREADs - $JOB_DIR/assembly_markdup.bam
samtools index $JOB_DIR/assembly_markdup.bam


