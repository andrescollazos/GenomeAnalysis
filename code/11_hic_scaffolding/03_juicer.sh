#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -J juicer
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/11_scaffolding/03_juicer/%x.%j.out

module load bioinfo-tools
module load samtools/1.20

THREADS=8

SRC_DIR="/home/anco6881/genome_analysis/GenomeAnalysis/analyses/11_scaffolding/02_scaffolding"
JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/11_scaffolding/03_juicer"

CONTIGS="$SRC_DIR/yahs.out_scaffolds_final.fa"
BIN_FILE="$SRC_DIR/yahs.out.bin"
SCAFFOLDS_FINAL="$SRC_DIR/yahs.out_scaffolds_final.agp"

# Step 1: Index the FASTA (needed for AGP parsing and size extraction)
cd $SRC_DIR
samtools faidx $CONTIGS

cd $JOB_DIR

# Step 2: Create sorted alignment file for juicer_tools input
/proj/uppmax2025-3-3/Genome_Analysis/yahs/juicer pre \
    $BIN_FILE $SCAFFOLDS_FINAL ${CONTIGS}.fai | \
    sort -k2,2d -k6,6d -T ./ --parallel=$THREADS -S32G | \
    awk 'NF' > alignments_sorted.txt

# Step 3: Extract scaffold sizes from the FASTA index
cut -f1,2 ${CONTIGS}.fai > scaffolds_final.chrom.sizes

# Step 4: Generate the .hic file for Juicebox visualization
java -Xmx32G -jar /proj/uppmax2025-3-3/Genome_Analysis/yahs/juicer_tools_1.22.01.jar pre \
    alignments_sorted.txt yahs_output.hic scaffolds_final.chrom.sizes

