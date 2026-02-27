#!/bin/bash
# Create a gene-level lookup file from a GENCODE GTF
# Usage: bash make_gene_lookup.sh <gencode_gtf> <output_file>
#
# Example:
#   bash make_gene_lookup.sh gencode.v49.annotation.gtf gencode_v49_genes.gtf
#   bash make_gene_lookup.sh gencode.v50.annotation.gtf gencode_v50_genes.gtf

set -e

if [[ $# -ne 2 ]]; then
    echo "Usage: bash make_gene_lookup.sh <gencode_gtf> <output_file>"
    echo "  <gencode_gtf>  Full GENCODE GTF (uncompressed or .gz)"
    echo "  <output_file>   Output gene-only GTF"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

echo "Extracting gene features from: $INPUT"

if [[ "$INPUT" == *.gz ]]; then
    zcat "$INPUT" | awk -F'\t' '$3 == "gene"' > "$OUTPUT"
else
    awk -F'\t' '$3 == "gene"' "$INPUT" > "$OUTPUT"
fi

NGENES=$(wc -l < "$OUTPUT")
SIZE=$(du -h "$OUTPUT" | cut -f1)
echo "Done: $NGENES genes -> $OUTPUT ($SIZE)"
