#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -J rna_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/06_rna_mapping/02_bam/%x.%j.out

# Load Modules
module load bioinfo-tools
module load samtools/1.20
module load bwa/0.7.18

export SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis/analyses/06_rna_mapping"
export JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/06_rna_mapping/02_bam"
export SAM_FILES_PATH="$SRC_DIR/01_mapping"

mkdir -p $JOB_DIR
cd "$JOB_DIR"

for sam_file in $SAM_FILES_PATH/*.sam
do
    base_name=$(basename "$sam_file" .sam)
    samtools view -bS "$sam_file" | samtools sort -o "${base_name}.sorted.bam"
    samtools index "${base_name}.sorted.bam"
done

ln -s "$JOB_DIR"/*.bam "$SRC_DIR"/02_bam