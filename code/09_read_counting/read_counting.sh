#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J read_counting
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/09_read_counting/%x.%j.out

module load bioinfo-tools subread/2.0.3

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/09_read_counting"

# Input Files
GTF_FILE="$SRC_DIR/analyses/07_structural_annotation/braker.gtf"
BAM_PATH="$SRC_DIR/analyses/06_rna_mapping/02_bam"
BAM_FILES=${BAM_PATH}/*.bam

cd "$JOB_DIR"

featureCounts -p --countReadPairs \
  -T 8 \
  -a "$GTF_FILE" -F GTF -t exon -g gene_id \
  -o "${JOB_DIR}/counts.txt" \
  $BAM_FILES
