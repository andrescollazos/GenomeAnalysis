#!/bin/bash

# Receive GTF file from command line
GTF_FILE="$1"

# Check if input was provided
if [[ -z "$GTF_FILE" ]]; then
    echo "Usage: $0 file.gtf"
    exit 1
fi

# Create output filename by removing .gtf and adding .statistics.txt
OUT_FILE="${GTF_FILE%.gtf}.statistics.txt"

# Create (or overwrite) the output file
{
echo "Statistics for $GTF_FILE"
echo "-------------------------------"

# Number of genes
genes=$(awk '$3 == "gene"' "$GTF_FILE" | wc -l)
echo "Number of genes: $genes"

# Number of transcripts (transcript features)
transcripts=$(awk '$3 == "transcript"' "$GTF_FILE" | wc -l)
echo "Number of transcripts: $transcripts"

# Number of exons
exons=$(awk '$3 == "exon"' "$GTF_FILE" | wc -l)
echo "Number of exons: $exons"

# Number of CDS
cds=$(awk '$3 == "CDS"' "$GTF_FILE" | wc -l)
echo "Number of CDS: $cds"

# Average exons per transcript
if [[ $transcripts -gt 0 ]]; then
    avg_exons=$(echo "scale=2; $exons / $transcripts" | bc)
    echo "Average exons per transcript: $avg_exons"
fi

# Average CDS per transcript
if [[ $transcripts -gt 0 ]]; then
    avg_cds=$(echo "scale=2; $cds / $transcripts" | bc)
    echo "Average CDS per transcript: $avg_cds"
fi
} > "$OUT_FILE"

echo "Statistics saved to $OUT_FILE"
