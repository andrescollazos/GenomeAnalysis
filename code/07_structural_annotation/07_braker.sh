#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 4:00:00
#SBATCH -J struc_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/07_structural_annotation/%x.%j.out

# Load Modules
module load bioinfo-tools braker/2.1.6

# Add Gm key
export HOME="/home/anco6881"
export GM_KEY="$HOME/.gm_key"

# Before running this script I must:
# module load bioinfo-tools braker/2.1.6
# Source Augustus config folder into the working directory (Check the script using cat $(echo $AUGUSTUS_CONFIG_COPY))
# cd $BACK_DIR/07_structural_annotation
# source $AUGUSTUS_CONFIG_COPY
# Provide permissions to the species subfolder: 
# chmod a+w -R augustus_config/species
# Copy the license key to $HOME
# cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/07_structural_annotation"

MASKED_GENOME="$SRC_DIR/analyses/05_masking/02_masking/polished_assembly.fasta.masked"
BAM_PATH="$SRC_DIR/analyses/06_rna_mapping/02_bam"
BAM_FILES=$(echo ${BAM_PATH}/*.bam | sed 's/ /,/g')

A_CONFIG_PATH="$JOB_DIR/augustus_config"

cd "$JOB_DIR"

braker.pl \
--genome="$MASKED_GENOME" \
--bam="$BAM_FILES" \
--softmasking=1 \
--species="niphotrichum_japonicum" \
--gff3 \
--cores=12 \
--AUGUSTUS_CONFIG_PATH="$A_CONFIG_PATH" \
--AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin \
--AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts \
--GENEMARK_PATH=/sw/bioinfo/GeneMark/4.68-es/snowy
