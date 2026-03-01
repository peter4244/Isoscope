#!/usr/bin/env bash
#
# extract_shortread_junctions.sh
#
# Extract junction-spanning reads from short-read BAMs for a genomic region and
# generate an IGV session. Filters reads to those containing splice junctions
# (N in CIGAR), producing small BAMs suitable for IGV visualization on hg38.
#
# Always produces both per-sample and merged BAMs in separate subdirectories:
#   <output-dir>/individual/   Per-sample junction BAMs
#   <output-dir>/merged/       Cell_type × treatment merged junction BAMs
#   <output-dir>/igv_session.xml   References merged BAMs by default
#   <output-dir>/gene_annotation.gtf   If --gtf provided
#
# Usage examples:
#   bash extract_shortread_junctions.sh --region chr18:8705271-8832780 \
#       --bam-dir /path/to/nfcore/results/ --output-dir ./sr_junctions/
#
#   bash extract_shortread_junctions.sh --gene MTCL1 \
#       --gencode-gtf /path/to/gencode.v47.annotation.gtf.gz \
#       --bam-dir /path/to/nfcore/results/ --bam-suffix .markdup.bam

set -euo pipefail

# ---- Defaults ----
REGION=""
GENE=""
GENCODE_GTF=""
BAM_DIR=""
BAM_FILES=()
BAM_SUFFIX=".markdup.bam"
OUTPUT_DIR=""
GTF=""
PAD=5000
DISPLAY_MODE="SQUISHED"
SAMPLE_REGEX='^Sample[0-9]+_(.+)_(DMSO|Smg1i)$'

# ---- Usage ----
usage() {
    cat <<'EOF'
Usage: extract_shortread_junctions.sh [OPTIONS]

Extract junction-spanning reads from short-read BAMs for a genomic region
and generate an IGV session. Always produces both individual and merged BAMs.

Region (one required):
  --region CHR:START-END   Genomic region (e.g., chr18:8705271-8832780)
  --gene GENE_NAME         Gene name to look up (requires --gencode-gtf)

BAM input (one required):
  --bam-dir DIR            Directory to search recursively for BAM files
  --bam FILE               Specific BAM file (repeatable)

Options:
  --bam-suffix SUFFIX      BAM filename suffix to match (default: .markdup.bam)
  --sample-regex REGEX     Regex for sample name parsing (default: ^Sample[0-9]+_(.+)_(DMSO|Smg1i)$)
  --output-dir DIR         Output directory (default: runs/<GENE>/ when --gene is used)
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
        --region)       REGION="$2"; shift 2 ;;
        --gene)         GENE="$2"; shift 2 ;;
        --gencode-gtf)  GENCODE_GTF="$2"; shift 2 ;;
        --bam-dir)      BAM_DIR="$2"; shift 2 ;;
        --bam)          BAM_FILES+=("$2"); shift 2 ;;
        --bam-suffix)   BAM_SUFFIX="$2"; shift 2 ;;
        --sample-regex) SAMPLE_REGEX="$2"; shift 2 ;;
        --output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
        --gtf)          GTF="$2"; shift 2 ;;
        --pad)          PAD="$2"; shift 2 ;;
        --display-mode) DISPLAY_MODE="$2"; shift 2 ;;
        -h|--help)      usage 0 ;;
        *)              echo "ERROR: Unknown option: $1"; usage 1 ;;
    esac
done

# ---- Validate dependencies ----
if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools is required but not found in PATH"
    echo "  Hint: conda activate /udd/repjc/.conda/envs/lr_igv"
    exit 1
fi

# ---- Detect samtools version for junction filtering method ----
SAMTOOLS_VERSION=$(samtools --version | head -1 | awk '{print $2}')
SAMTOOLS_MAJOR=$(echo "$SAMTOOLS_VERSION" | cut -d. -f1)
SAMTOOLS_MINOR=$(echo "$SAMTOOLS_VERSION" | cut -d. -f2)

USE_EXPR_FILTER=false
if [[ "$SAMTOOLS_MAJOR" -gt 1 ]] || [[ "$SAMTOOLS_MAJOR" -eq 1 && "$SAMTOOLS_MINOR" -ge 12 ]]; then
    USE_EXPR_FILTER=true
fi
echo "samtools ${SAMTOOLS_VERSION} detected (expression filter: ${USE_EXPR_FILTER})"

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
    if [[ -n "$GENE" ]]; then
        OUTPUT_DIR="runs/${GENE}"
        echo "Defaulting output directory to: ${OUTPUT_DIR}/"
    else
        echo "ERROR: --output-dir is required (or use --gene to auto-default to runs/<GENE>/)"
        usage 1
    fi
fi

# ---- Step 1: Resolve gene coordinates ----
if [[ -n "$GENE" ]]; then
    echo "Looking up coordinates for gene: $GENE"

    if [[ ! -f "$GENCODE_GTF" ]]; then
        echo "ERROR: GENCODE GTF not found: $GENCODE_GTF"
        exit 1
    fi

    # Determine if file is compressed
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
REGION_CHR=$(echo "$REGION" | cut -d: -f1)
REGION_COORDS=$(echo "$REGION" | cut -d: -f2)
REGION_START=$(echo "$REGION_COORDS" | cut -d- -f1)
REGION_END=$(echo "$REGION_COORDS" | cut -d- -f2)

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
    # Recursive search to handle nfcore-rnaseq subdirectory layout
    while IFS= read -r -d '' bam; do
        BAM_FILES+=("$bam")
    done < <(find "$BAM_DIR" -name "*${BAM_SUFFIX}" -print0 | sort -z)
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found (suffix: ${BAM_SUFFIX})"
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM file(s)"

# ---- Step 3: Check chromosome naming convention ----
FIRST_BAM="${BAM_FILES[0]}"
BAM_HAS_CHR=$(samtools idxstats "$FIRST_BAM" 2>/dev/null \
    | awk '$1 ~ /^chr[0-9]/' | head -1)

REGION_HAS_CHR=false
if [[ "$REGION_CHR" == chr* ]]; then
    REGION_HAS_CHR=true
fi

if [[ "$REGION_HAS_CHR" == true && -z "$BAM_HAS_CHR" ]]; then
    QUERY_CHR="${REGION_CHR#chr}"
    echo "WARNING: BAM uses non-chr naming. Stripping 'chr' prefix for queries."
elif [[ "$REGION_HAS_CHR" == false && -n "$BAM_HAS_CHR" ]]; then
    QUERY_CHR="chr${REGION_CHR}"
    echo "WARNING: BAM uses chr-prefixed naming. Adding 'chr' prefix for queries."
else
    QUERY_CHR="$REGION_CHR"
fi

QUERY_REGION="${QUERY_CHR}:${PADDED_START}-${PADDED_END}"

# ---- Step 4: Create output directories and extract junction reads ----
INDIVIDUAL_DIR="${OUTPUT_DIR}/individual"
MERGED_DIR="${OUTPUT_DIR}/merged"
mkdir -p "$INDIVIDUAL_DIR" "$MERGED_DIR"

# Strip BAM_SUFFIX to get the stem for sample name parsing
# e.g., ".markdup.bam" → strip both ".markdup" and ".bam"
STEM_SUFFIX="${BAM_SUFFIX%.bam}"  # e.g., ".markdup"

EXTRACTED_BAMS=()
EMPTY_BAMS=()

for bam in "${BAM_FILES[@]}"; do
    bam_basename=$(basename "$bam")

    # Validate BAM index exists
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "  WARNING: No index for ${bam_basename}, attempting to create..."
        samtools index "$bam"
    fi

    # Derive output name: strip original suffix, add .junction.bam
    out_stem="${bam_basename%${BAM_SUFFIX}}"
    out_name="${out_stem}.junction.bam"
    out_bam="${INDIVIDUAL_DIR}/${out_name}"

    echo "  Extracting junctions: ${bam_basename}"

    if [[ "$USE_EXPR_FILTER" == true ]]; then
        # samtools >= 1.12: use expression filter for junction reads
        samtools view -b -h -e 'cigar =~ "N"' "$bam" "$QUERY_REGION" > "$out_bam"
    else
        # Fallback: awk-based CIGAR filtering
        samtools view -h "$bam" "$QUERY_REGION" \
            | awk '$1 ~ /^@/ || $6 ~ /[0-9]+N/' \
            | samtools view -bS - > "$out_bam"
    fi

    # Check if extraction produced reads
    READ_COUNT=$(samtools view -c "$out_bam")
    if [[ "$READ_COUNT" -eq 0 ]]; then
        echo "    WARNING: 0 junction reads from ${bam_basename}"
        EMPTY_BAMS+=("$out_name")
    else
        echo "    ${READ_COUNT} junction reads"
    fi

    samtools index "$out_bam"
    EXTRACTED_BAMS+=("$out_name")
done

echo ""
echo "Extracted ${#EXTRACTED_BAMS[@]} junction BAM(s) to ${INDIVIDUAL_DIR}/"
if [[ ${#EMPTY_BAMS[@]} -gt 0 ]]; then
    echo "WARNING: ${#EMPTY_BAMS[@]} BAM(s) had 0 junction reads: ${EMPTY_BAMS[*]}"
fi

# ---- Step 4b: Merge BAMs by group (cell_type × treatment) ----
echo ""
echo "Merging BAMs by (cell_type, treatment) group..."

GROUP_LINES=""
for bam_name in "${EXTRACTED_BAMS[@]}"; do
    # Strip .junction.bam to get sample stem
    sample_id="${bam_name%.junction.bam}"

    # Also strip any remaining suffix components (e.g., ".markdup" already gone)
    ct=""
    treatment=""
    if [[ "$sample_id" =~ $SAMPLE_REGEX ]]; then
        middle="${BASH_REMATCH[1]}"
        treatment="${BASH_REMATCH[2]}"
        donor="${middle##*_}"
        ct="${middle%_${donor}}"
    fi
    if [[ -n "$ct" && -n "$treatment" ]]; then
        GROUP_LINES="${GROUP_LINES}${ct}|${treatment}|${bam_name}
"
    else
        echo "  WARNING: Could not parse sample name '${bam_name}', skipping merge for this file"
    fi
done

SORTED_GROUPS=$(echo "$GROUP_LINES" | sort -t'|' -k1,1 -k2,2)
UNIQUE_GROUPS=$(echo "$SORTED_GROUPS" | cut -d'|' -f1,2 | sort -u)

MERGED_BAMS=()
MERGED_DISPLAY_LINES=""

while IFS='|' read -r grp_ct grp_treat; do
    [[ -z "$grp_ct" ]] && continue

    GRP_BAM_LIST=()
    while IFS='|' read -r line_ct line_treat line_bam; do
        if [[ "$line_ct" == "$grp_ct" && "$line_treat" == "$grp_treat" ]]; then
            GRP_BAM_LIST+=("$line_bam")
        fi
    done <<< "$SORTED_GROUPS"

    N=${#GRP_BAM_LIST[@]}
    MERGED_NAME="${grp_ct}_${grp_treat}.junction.bam"
    DISPLAY_NAME="${grp_ct} ${grp_treat} (n=${N})"

    if [[ $N -eq 1 ]]; then
        cp "${INDIVIDUAL_DIR}/${GRP_BAM_LIST[0]}" "${MERGED_DIR}/${MERGED_NAME}"
        cp "${INDIVIDUAL_DIR}/${GRP_BAM_LIST[0]}.bai" "${MERGED_DIR}/${MERGED_NAME}.bai"
        echo "  ${DISPLAY_NAME}: copied ${GRP_BAM_LIST[0]}"
    else
        MERGE_INPUTS=()
        for b in "${GRP_BAM_LIST[@]}"; do
            MERGE_INPUTS+=("${INDIVIDUAL_DIR}/${b}")
        done
        samtools merge -c -f "${MERGED_DIR}/${MERGED_NAME}" "${MERGE_INPUTS[@]}"
        samtools index "${MERGED_DIR}/${MERGED_NAME}"
        echo "  ${DISPLAY_NAME}: merged ${N} BAMs"
    fi

    MERGED_BAMS+=("$MERGED_NAME")
    MERGED_DISPLAY_LINES="${MERGED_DISPLAY_LINES}${grp_ct}|${grp_treat}|${MERGED_NAME}|${DISPLAY_NAME}
"
done <<< "$UNIQUE_GROUPS"

echo "Merged into ${#MERGED_BAMS[@]} group(s) in ${MERGED_DIR}/"

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
IGV_LOCUS="${REGION_CHR}:${PADDED_START}-${PADDED_END}"
SESSION_FILE="${OUTPUT_DIR}/igv_session.xml"

SORTED_BAMS=$(echo "$MERGED_DISPLAY_LINES" | while IFS='|' read -r grp_ct grp_treat bam_name display_name; do
    [[ -z "$grp_ct" ]] && continue
    treat_order=2
    if [[ "$grp_treat" == "DMSO" ]]; then treat_order=0; fi
    if [[ "$grp_treat" == "Smg1i" ]]; then treat_order=1; fi
    echo "${grp_ct}||${treat_order}|merged/${bam_name}|${display_name}"
done | sort -t'|' -k1,1 -k3,3n)

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
echo "  Tracks: ${#MERGED_BAMS[@]} merged group(s)"
echo "  Individual BAMs: ${INDIVIDUAL_DIR}/ (${#EXTRACTED_BAMS[@]} files)"
echo "  Merged BAMs: ${MERGED_DIR}/ (${#MERGED_BAMS[@]} files)"
if [[ -n "$GTF_FILENAME" ]]; then
    echo "  Annotation: ${GTF_FILENAME}"
fi
echo ""
echo "To open: File > Open Session in IGV, or:"
echo "  igv ${SESSION_FILE}"
