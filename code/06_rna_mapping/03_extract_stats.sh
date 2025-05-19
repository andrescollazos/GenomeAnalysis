#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -J bam_coverage
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andres-felipe.collazos-rozo.6881@student.uu.se
#SBATCH --output=/home/anco6881/genome_analysis/GenomeAnalysis/analyses/06_rna_mapping/02_bam/%x.%j.out

module load bioinfo-tools
module load samtools/1.20

JOB_DIR="/proj/uppmax2025-3-3/nobackup/work/anco6881/06_rna_mapping/02_bam"
cd $JOB_DIR

echo "Processing Control_1"
samtools coverage -m Control_1.sorted.bam -o Control_1.coverage.txt
echo "Processing Control_2"
samtools coverage -m Control_2.sorted.bam -o Control_2.coverage.txt
echo "Processing Control_3"
samtools coverage -m Control_3.sorted.bam -o Control_3.coverage.txt
echo "Processing Heat_treated_42_12h_1"
samtools coverage -m Heat_treated_42_12h_1.sorted.bam -o Heat_treated_42_12h_1.coverage.txt
echo "Processing Heat_treated_42_12h_2"
samtools coverage -m Heat_treated_42_12h_2.sorted.bam -o Heat_treated_42_12h_2.coverage.txt
echo "Processing Heat_treated_42_12h_3"
samtools coverage -m Heat_treated_42_12h_3.sorted.bam -o Heat_treated_42_12h_3.coverage.txt

for file in *.coverage.txt
do
    echo -n "$file: "; grep "Mean coverage:" "$file" | \
    sed -E 's/.*Mean coverage:[[:space:]]+([0-9.]+)x/\1/' | \
    awk '{sum+=$1} END {print sum/NR}'
done
