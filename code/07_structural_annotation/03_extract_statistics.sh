# Load Modules
#module load bioinfo-tools AGAT/1.3.2

JOB_DIR="/home/anco6881/genome_analysis/GenomeAnalysis/analyses/07_structural_annotation/"
GFF="$JOB_DIR/braker_standardized.gff3"

cd $JOB_DIR

# Convert GFF file to .tsv
agat_convert_sp_gff2tsv.pl --gff $GFF -o features.tsv

# --------------------------------------------------------------------
# Extract Gene lengths
awk -F'\t' '$3 == "gene" {print $5 - $4 + 1}' features.tsv > gene_lengths.txt

# --------------------------------------------------------------------
# Extract counts per exon
awk -F'\t' '$3 == "exon" {print $10}' features.tsv | sort | uniq -c > exon_counts_per_mrna.txt
awk '{print $1}' exon_counts_per_mrna.txt > exon_counts_only.txt

