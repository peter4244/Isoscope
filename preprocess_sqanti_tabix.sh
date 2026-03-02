#!/bin/bash
# Preprocess GENCODE and SQANTI GTF files with tabix indexing
# This is a ONE-TIME preprocessing step (~10-15 minutes)
# Creates compressed + indexed GTFs for <1 second queries by genomic coordinates

set -e  # Exit on error

# Load file paths from config (copy config.example.sh -> config.sh and edit)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/config.sh" ]]; then
    echo "ERROR: config.sh not found. Copy config.example.sh to config.sh and edit paths."
    exit 1
fi
source "${SCRIPT_DIR}/config.sh"

# Derived paths (from GENCODE_GTF and SQANTI_GTF)
GENCODE_SORTED="${GENCODE_GTF%.gtf}.sorted.gtf"
GENCODE_GZ="${GENCODE_GTF%.gtf}.sorted.gtf.gz"

SQANTI_SORTED="${SQANTI_GTF%.gtf}.sorted.gtf"
SQANTI_GZ="${SQANTI_GTF%.gtf}.sorted.gtf.gz"

# Check if input files exist
if [ ! -f "$GENCODE_GTF" ]; then
    echo "ERROR: GENCODE GTF not found: $GENCODE_GTF"
    exit 1
fi

if [ ! -f "$SQANTI_GTF" ]; then
    echo "ERROR: SQANTI GTF not found: $SQANTI_GTF"
    exit 1
fi

# Check if tools are available
if ! command -v bgzip &> /dev/null; then
    echo "ERROR: bgzip not found. Install with:"
    echo "  conda install -c bioconda htslib"
    echo "  OR"
    echo "  brew install htslib"
    exit 1
fi

if ! command -v tabix &> /dev/null; then
    echo "ERROR: tabix not found. Install with:"
    echo "  conda install -c bioconda htslib"
    echo "  OR"
    echo "  brew install htslib"
    exit 1
fi

echo "======================================================================"
echo "GTF Tabix Preprocessing for Gene Isoform Annotation Pipeline"
echo "======================================================================"
echo "This will create tabix-indexed versions of:"
echo "  1. GENCODE GTF (2.9 GB)"
echo "  2. SQANTI GTF (1.0 GB)"
echo ""
echo "Estimated time: 10-15 minutes"
echo "======================================================================"
echo ""

# Function to process a GTF file
process_gtf() {
    local INPUT=$1
    local SORTED=$2
    local OUTPUT=$3
    local NAME=$4

    echo "Processing $NAME..."
    echo "  Input:  $INPUT"
    echo "  Output: $OUTPUT"
    echo ""

    # Check if output already exists
    if [ -f "$OUTPUT" ] && [ -f "${OUTPUT}.tbi" ]; then
        echo "  WARNING: Output files already exist. Skipping."
        echo "  (Delete $OUTPUT to reprocess)"
        echo ""
        return 0
    fi

    # Step 1: Sort
    echo "  [1/3] Sorting by genomic coordinates..."
    grep -v "^#" "$INPUT" | sort -k1,1 -k4,4n > "$SORTED"

    # Step 2: Compress
    echo "  [2/3] Compressing with bgzip..."
    bgzip -c "$SORTED" > "$OUTPUT"

    # Step 3: Index
    echo "  [3/3] Creating tabix index..."
    tabix -p gff "$OUTPUT"

    # Clean up
    rm "$SORTED"

    echo "  ✓ Complete!"
    echo "    Compressed: $(du -h "$OUTPUT" | cut -f1)"
    echo "    Index:      $(du -h "${OUTPUT}.tbi" | cut -f1)"
    echo ""
}

# Process GENCODE GTF
process_gtf "$GENCODE_GTF" "$GENCODE_SORTED" "$GENCODE_GZ" "GENCODE v49"

# Process SQANTI GTF
process_gtf "$SQANTI_GTF" "$SQANTI_SORTED" "$SQANTI_GZ" "SQANTI PacBio"

echo "======================================================================"
echo "SUCCESS! All GTF files have been preprocessed."
echo "======================================================================"
echo ""
echo "Output files:"
ls -lh "$GENCODE_GZ" "${GENCODE_GZ}.tbi"
ls -lh "$SQANTI_GZ" "${SQANTI_GZ}.tbi"
echo ""
echo "The gene_isoform_annotation.R script will now use these indexed files"
echo "for fast (<1 second) queries."
echo ""
