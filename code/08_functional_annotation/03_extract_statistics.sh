JOB_DIR="../../analyses/08_functional_annotation/02_eggNOG_mapper"
ANNOTATION="$JOB_DIR/n_japonicum.emapper.annotations"

# Total predicted proteins (non-comment lines)
echo "Total predicted proteins:" && grep -v '^#' $ANNOTATION | wc -l

# With COG category (column 7)
echo "Proteins with assigned COG category:" && awk -F'\t' '!/^#/ && $7 != "" && $7 != "-"' $ANNOTATION | wc -l

# With GO terms (column 10)
echo "Proteins with assigned GO terms:" && awk -F'\t' '!/^#/ && $10 != "" && $10 != "-"' $ANNOTATION | wc -l

# With EC numbers (column 11)
echo "Proteins with assigned EC numbers:" && awk -F'\t' '!/^#/ && $11 != "" && $11 != "-"' $ANNOTATION | wc -l

# With KEGG orthologs (column 12)
echo "Proteins with assigned KEGG orthologs:" && awk -F'\t' '!/^#/ && $12 != "" && $12 != "-"' $ANNOTATION | wc -l

# Extract COG category column (column 7), flatten multi-letter entries
cut -f7 $ANNOTATION | grep -v '^#' | grep -v '^\s*$' | tr ',' '\n' | tr -d ' ' | sort | uniq -c | sort -nr > "$JOB_DIR/cog_counts.txt"
