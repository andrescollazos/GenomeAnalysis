#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J scaffolding
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/11_scaffolding/02_scaffolding/%x.%j.out

module load bioinfo-tools

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/11_scaffolding/02_scaffolding"

REF_GENOME="$SRC_DIR/data/polished_assembly/polished_assembly.fasta"
HIC_TO_CONTIGS="$SRC_DIR/analyses/11_scaffolding/01_hic_mapping/hic_to_contigs.bam"

cd $JOB_DIR
/proj/uppmax2025-3-3/Genome_Analysis/yahs/yahs $REF_GENOME $HIC_TO_CONTIGS