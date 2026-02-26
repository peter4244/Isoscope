#!/usr/bin/env bash
#
# extract_gene_region.sh
#
# Extract a genomic region from genome-aligned BAMs and generate an IGV session.
# Produces small, self-contained BAM slices suitable for IGV visualization on hg38.
#
# Usage examples:
#   bash extract_gene_region.sh --region chr18:8705271-8832780 \
#       --bam-dir /path/to/bams/ --output-dir ./igv_genome/
#
#   bash extract_gene_region.sh --region chr18:8705271-8832780 \
#       --bam /path/to/sample1.bam --bam /path/to/sample2.bam \
#       --output-dir ./igv_genome/ --gtf annotation.gtf
#
#   bash extract_gene_region.sh --gene MTCL1 \
#       --gencode-gtf /path/to/gencode.v47.annotation.gtf.gz \
#       --bam-dir /path/to/bams/ --output-dir ./igv_genome/

set -euo pipefail

# ---- Defaults ----
REGION=""
GENE=""
GENCODE_GTF=""
BAM_DIR=""
BAM_FILES=()
OUTPUT_DIR=""
GTF=""
PAD=5000
DISPLAY_MODE="SQUISHED"

# ---- Usage ----
usage() {
    cat <<'EOF'
Usage: extract_gene_region.sh [OPTIONS]

Extract a genomic region from genome-aligned BAMs and generate an IGV session.

Region (one required):
  --region CHR:START-END   Genomic region (e.g., chr18:8705271-8832780)
  --gene GENE_NAME         Gene name to look up (requires --gencode-gtf)

BAM input (one required):
  --bam-dir DIR            Directory containing .bam files
  --bam FILE               Specific BAM file (repeatable)

Options:
  --output-dir DIR         Output directory (required)
  --gtf FILE               Gene annotation GTF to include in session
  --gencode-gtf FILE       Tabix-indexed or bgzipped GENCODE GTF for gene lookup
  --pad N                  Padding around region in bp (default: 5000)
  --display-mode MODE      IGV display mode: SQUISHED, EXPANDED, COLLAPSED (default: SQUISHED)
  -h, --help               Show this help
EOF
    exit "${1:-0}"
}

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --region)      REGION="$2"; shift 2 ;;
        --gene)        GENE="$2"; shift 2 ;;
        --gencode-gtf) GENCODE_GTF="$2"; shift 2 ;;
        --bam-dir)     BAM_DIR="$2"; shift 2 ;;
        --bam)         BAM_FILES+=("$2"); shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
        --gtf)         GTF="$2"; shift 2 ;;
        --pad)         PAD="$2"; shift 2 ;;
        --display-mode) DISPLAY_MODE="$2"; shift 2 ;;
        -h|--help)     usage 0 ;;
        *)             echo "ERROR: Unknown option: $1"; usage 1 ;;
    esac
done

# ---- Validate dependencies ----
if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools is required but not found in PATH"
    exit 1
fi

# ---- Validate required arguments ----
if [[ -z "$REGION" && -z "$GENE" ]]; then
    echo "ERROR: Either --region or --gene is required"
    usage 1
fi

if [[ -n "$GENE" && -z "$GENCODE_GTF" ]]; then
    echo "ERROR: --gene requires --gencode-gtf"
    exit 1
fi

if [[ -z "$BAM_DIR" && ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: Either --bam-dir or --bam is required"
    usage 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --output-dir is required"
    usage 1
fi

# ---- Step 1: Resolve gene coordinates ----
if [[ -n "$GENE" ]]; then
    echo "Looking up coordinates for gene: $GENE"

    if [[ ! -f "$GENCODE_GTF" ]]; then
        echo "ERROR: GENCODE GTF not found: $GENCODE_GTF"
        exit 1
    fi

    # Determine if file is compressed
    # Note: grep -P is not available on macOS; use awk for tab-delimited field matching
    if [[ "$GENCODE_GTF" == *.gz ]]; then
        GENE_LINE=$(zcat "$GENCODE_GTF" \
            | awk -F'\t' '$3 == "gene"' \
            | grep "gene_name \"${GENE}\"" \
            | head -1)
    else
        GENE_LINE=$(awk -F'\t' '$3 == "gene"' "$GENCODE_GTF" \
            | grep "gene_name \"${GENE}\"" \
            | head -1)
    fi

    if [[ -z "$GENE_LINE" ]]; then
        echo "ERROR: Gene '$GENE' not found in $GENCODE_GTF"
        exit 1
    fi

    GENE_CHR=$(echo "$GENE_LINE" | awk -F'\t' '{print $1}')
    GENE_START=$(echo "$GENE_LINE" | awk -F'\t' '{print $4}')
    GENE_END=$(echo "$GENE_LINE" | awk -F'\t' '{print $5}')

    echo "  Found: ${GENE_CHR}:${GENE_START}-${GENE_END}"
    REGION="${GENE_CHR}:${GENE_START}-${GENE_END}"
fi

# ---- Apply padding ----
# Parse region components
REGION_CHR=$(echo "$REGION" | cut -d: -f1)
REGION_COORDS=$(echo "$REGION" | cut -d: -f2)
REGION_START=$(echo "$REGION_COORDS" | cut -d- -f1)
REGION_END=$(echo "$REGION_COORDS" | cut -d- -f2)

# Apply padding (floor at 1)
PADDED_START=$(( REGION_START - PAD ))
if (( PADDED_START < 1 )); then
    PADDED_START=1
fi
PADDED_END=$(( REGION_END + PAD ))

PADDED_REGION="${REGION_CHR}:${PADDED_START}-${PADDED_END}"
echo "Region (with ${PAD}bp padding): ${PADDED_REGION}"

# ---- Step 2: Discover BAM files ----
if [[ -n "$BAM_DIR" ]]; then
    if [[ ! -d "$BAM_DIR" ]]; then
        echo "ERROR: BAM directory not found: $BAM_DIR"
        exit 1
    fi
    while IFS= read -r -d '' bam; do
        BAM_FILES+=("$bam")
    done < <(find "$BAM_DIR" -maxdepth 1 -name '*.bam' -print0 | sort -z)
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found"
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM file(s)"

# ---- Step 3: Check chromosome naming convention ----
# Compare BAM SQ headers with the region's chromosome name
FIRST_BAM="${BAM_FILES[0]}"
BAM_HAS_CHR=$(samtools idxstats "$FIRST_BAM" 2>/dev/null \
    | awk '$1 ~ /^chr[0-9]/' | head -1)

REGION_HAS_CHR=false
if [[ "$REGION_CHR" == chr* ]]; then
    REGION_HAS_CHR=true
fi

if [[ "$REGION_HAS_CHR" == true && -z "$BAM_HAS_CHR" ]]; then
    # Region has chr prefix but BAM does not — strip it
    QUERY_CHR="${REGION_CHR#chr}"
    echo "WARNING: BAM uses non-chr naming. Stripping 'chr' prefix for queries."
elif [[ "$REGION_HAS_CHR" == false && -n "$BAM_HAS_CHR" ]]; then
    # Region lacks chr prefix but BAM has it — add it
    QUERY_CHR="chr${REGION_CHR}"
    echo "WARNING: BAM uses chr-prefixed naming. Adding 'chr' prefix for queries."
else
    QUERY_CHR="$REGION_CHR"
fi

QUERY_REGION="${QUERY_CHR}:${PADDED_START}-${PADDED_END}"

# ---- Step 4: Create output directory and extract regions ----
mkdir -p "$OUTPUT_DIR"

EXTRACTED_BAMS=()
EMPTY_BAMS=()

for bam in "${BAM_FILES[@]}"; do
    bam_basename=$(basename "$bam")

    # Validate BAM index exists
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "  WARNING: No index for ${bam_basename}, attempting to create..."
        samtools index "$bam"
    fi

    out_bam="${OUTPUT_DIR}/${bam_basename}"
    echo "  Extracting: ${bam_basename}"
    samtools view -b -h "$bam" "$QUERY_REGION" > "$out_bam"

    # Check if extraction produced reads
    READ_COUNT=$(samtools view -c "$out_bam")
    if [[ "$READ_COUNT" -eq 0 ]]; then
        echo "    WARNING: 0 reads extracted from ${bam_basename} (region may not be covered)"
        EMPTY_BAMS+=("$bam_basename")
    fi

    samtools index "$out_bam"
    EXTRACTED_BAMS+=("$bam_basename")
done

echo ""
echo "Extracted ${#EXTRACTED_BAMS[@]} BAM(s) to ${OUTPUT_DIR}/"
if [[ ${#EMPTY_BAMS[@]} -gt 0 ]]; then
    echo "WARNING: ${#EMPTY_BAMS[@]} BAM(s) had 0 reads in the region: ${EMPTY_BAMS[*]}"
fi

# ---- Step 5: Copy GTF annotation if provided ----
GTF_FILENAME=""
if [[ -n "$GTF" ]]; then
    if [[ ! -f "$GTF" ]]; then
        echo "WARNING: GTF not found: $GTF (skipping)"
    else
        GTF_FILENAME="gene_annotation.gtf"
        cp "$GTF" "${OUTPUT_DIR}/${GTF_FILENAME}"
        echo "Copied annotation GTF to ${OUTPUT_DIR}/${GTF_FILENAME}"
    fi
fi

# ---- Step 6: Generate IGV session XML ----
#
# Track ordering: sorted by cell type, donor, then treatment (DMSO before Smg1i).
# Sample names like "Sample13_DD_017Q_DMSO" are parsed to extract a display name.
# For the IGV locus, use the original (chr-prefixed) region for hg38 compatibility.

IGV_LOCUS="${REGION_CHR}:${PADDED_START}-${PADDED_END}"
SESSION_FILE="${OUTPUT_DIR}/igv_session.xml"

# Build sorted BAM list with sort keys
# Parse sample name: Sample{N}_{CT}_{DONOR}_{TREATMENT}.aligned.bam
# Sort by: cell type, donor, treatment (DMSO < Smg1i)
SORTED_BAMS=$(for bam_name in "${EXTRACTED_BAMS[@]}"; do
    # Strip extensions to get sample ID
    sample_id="${bam_name%.aligned.bam}"
    sample_id="${sample_id%.bam}"

    # Extract components (Sample{N}_{CT}_{DONOR}_{TREATMENT})
    # Handle names like Sample1_DD_ALI_001V_Smg1i (cell type with underscore)
    ct=""
    donor=""
    treatment=""
    display_name="$sample_id"

    # Try to parse the sample name
    if [[ "$sample_id" =~ ^Sample[0-9]+_(.+)_(DMSO|Smg1i)$ ]]; then
        middle="${BASH_REMATCH[1]}"
        treatment="${BASH_REMATCH[2]}"
        # middle is like "DD_017Q" or "DD_ALI_001V"
        # donor is the last part (alphanumeric ID like 017Q, 001V)
        donor="${middle##*_}"
        ct="${middle%_${donor}}"
        display_name="${ct}_${donor}_${treatment}"
    fi

    # Sort key: ct, donor, treatment (DMSO=0, Smg1i=1)
    treat_order=2
    if [[ "$treatment" == "DMSO" ]]; then treat_order=0; fi
    if [[ "$treatment" == "Smg1i" ]]; then treat_order=1; fi

    echo "${ct}|${donor}|${treat_order}|${bam_name}|${display_name}"
done | sort -t'|' -k1,1 -k2,2 -k3,3n)

# Write the XML
{
    echo '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
    echo "<Session genome=\"hg38\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"${IGV_LOCUS}\" version=\"8\">"
    echo '    <Resources>'

    while IFS='|' read -r _ _ _ bam_name display_name; do
        echo "        <Resource path=\"${bam_name}\" name=\"${display_name}\"/>"
    done <<< "$SORTED_BAMS"

    if [[ -n "$GTF_FILENAME" ]]; then
        echo "        <Resource path=\"${GTF_FILENAME}\" name=\"Gene Annotation\"/>"
    fi

    echo '    </Resources>'
    echo '    <Panel name="DataPanel" height="600">'

    while IFS='|' read -r _ _ _ bam_name display_name; do
        cat <<TRACK
        <Track id="${bam_name}" name="${display_name}" displayMode="${DISPLAY_MODE}" showAllBases="true" visible="true">
            <RenderOptions/>
        </Track>
TRACK
    done <<< "$SORTED_BAMS"

    echo '    </Panel>'

    if [[ -n "$GTF_FILENAME" ]]; then
        echo '    <Panel name="FeaturePanel" height="120">'
        echo "        <Track id=\"${GTF_FILENAME}\" name=\"Gene Annotation\" displayMode=\"EXPANDED\" visible=\"true\">"
        echo '            <RenderOptions/>'
        echo '        </Track>'
        echo '    </Panel>'
        echo '    <PanelLayout dividerFractions="0.85"/>'
    else
        echo '    <PanelLayout dividerFractions="0.95"/>'
    fi

    echo '</Session>'
} > "$SESSION_FILE"

echo ""
echo "IGV session: ${SESSION_FILE}"
echo "  Genome: hg38"
echo "  Locus: ${IGV_LOCUS}"
echo "  Tracks: ${#EXTRACTED_BAMS[@]} BAM(s)"
if [[ -n "$GTF_FILENAME" ]]; then
    echo "  Annotation: ${GTF_FILENAME}"
fi
echo ""
echo "To open: File > Open Session in IGV, or:"
echo "  igv ${SESSION_FILE}"
