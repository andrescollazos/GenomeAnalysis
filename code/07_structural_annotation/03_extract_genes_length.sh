# Load Modules
module load bioinfo-tools AGAT/1.3.2

JOB_DIR="/home/anco6881/genome_analysis/GenomeAnalysis/analyses/07_structural_annotation/"
GFF="$JOB_DIR/braker_standardized.gff3"

agat_convert_sp_gff2tsv.pl --gff $GFF -o "$JOB_DIR/features.tsv"
awk -F'\t' '$3 == "gene" {print $5 - $4 + 1}' "$JOB_DIR/features.tsv" > "$JOB_DIR/gene_lengths.txt"
