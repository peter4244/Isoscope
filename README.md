# lr_igv_from_bam

A toolkit for PacBio long-read isoform analysis, providing:

- **Isoform Annotation Pipeline** — Extract comprehensive transcript-level annotations for any gene, combining GENCODE reference annotations with SQANTI PacBio long-read isoforms. Produces annotated TSVs, GTFs, FASTAs, and protein sequences with optional expression quantification and UniProt comparison.
- **BAM Visualization Tools** — Extract gene regions from genome-aligned (hg38) PacBio BAMs into small, self-contained files with auto-generated IGV sessions and sashimi plot support.

Developed for the NMD (Nonsense-Mediated Decay) lung cell line study but works with any GENCODE + SQANTI + PacBio dataset.

---

## Table of Contents

- [Setup](#setup)
  - [Prerequisites](#prerequisites)
  - [Configuration](#configuration)
  - [One-Time Preprocessing](#one-time-preprocessing)
  - [Required Files](#required-files)
- [Isoform Annotation Pipeline](#isoform-annotation-pipeline)
  - [Quick Start](#quick-start)
  - [Command-Line Options](#command-line-options)
  - [Expression Configuration](#expression-configuration)
  - [Output Files](#output-files)
  - [Sample Output](#sample-output)
- [BAM Visualization & IGV Sessions](#bam-visualization--igv-sessions)
  - [Quick Start](#quick-start-1)
  - [extract_gene_region.sh](#extract_gene_regionsh)
  - [Sashimi Plots with ggsashimi](#sashimi-plots-with-ggsashimi)
- [Reference](#reference)
  - [How the Annotation Pipeline Works](#how-the-annotation-pipeline-works)
  - [Protein Analysis Algorithms](#protein-analysis-algorithms)
  - [Data Flow Diagram](#data-flow-diagram)
  - [NMD Study Terminology](#nmd-study-terminology)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Version History](#version-history)
- [Support](#support)
- [License](#license)

---

## Setup

### Prerequisites

#### System Requirements

- **OS**: Linux or macOS
- **RAM**: 8 GB minimum, 16 GB recommended
- **Disk**: ~5 GB for reference files

#### Option A: Conda Environment (Recommended)

A self-contained conda environment bundles R, all R packages, and command-line tools so nothing needs to be compiled from source. This avoids compiler version issues common on HPC clusters.

```bash
# Set strict channel priority (avoids ABI conflicts between conda-forge and bioconda)
conda config --set channel_priority strict

# Create and activate the environment
conda env create -f environment.yml
conda activate lr_igv
```

After activation, all R packages, `tabix`, `bgzip`, and `samtools` are available — skip to [Configuration](#configuration).

#### Option B: Manual Installation

<details>
<summary>Click to expand manual installation steps</summary>

##### R (version >= 4.0)
Install from https://cran.r-project.org/

##### R Packages
```r
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "rtracklayer",      # GTF/GFF parsing
  "GenomicRanges",    # Genomic feature handling
  "edgeR",            # DGEList handling, CPM calculation
  "Biostrings"        # FASTA handling and translation (needed for --fasta and --protein)
))

# Install CRAN packages
install.packages(c(
  "dplyr",            # Data manipulation
  "tidyr",            # Data tidying
  "readr"             # Fast file reading
))
```

##### tabix and bgzip (part of htslib)

**Required for preprocessing and running the annotation pipeline.**

**Via Conda:**
```bash
conda install -c bioconda htslib
```

**Via Homebrew** (macOS):
```bash
brew install htslib
```

**Via apt** (Ubuntu/Debian):
```bash
sudo apt-get install tabix
```

**Purpose:** Enables fast coordinate-based GTF queries (<1 second) for both GENCODE and SQANTI files.

**Verification:**
```bash
which tabix && which bgzip && echo "Tools installed!" || echo "Need to install htslib"
```

##### samtools

**Required for BAM visualization tools** (`extract_gene_region.sh`).

```bash
conda install -c bioconda samtools   # or: brew install samtools
```

</details>

### Configuration

```bash
# 1. Clone the repo
git clone https://changit.bwh.harvard.edu/repjc/lr_igv_from_bam.git && cd lr_igv_from_bam

# 2. Activate the conda environment (if using Option A)
conda activate lr_igv

# 3. Create config files from templates and edit paths for your system
cp config.example.R config.R    # edit paths — used by gene_isoform_annotation.R
cp config.example.sh config.sh  # edit paths — used by preprocess_sqanti_tabix.sh
```

### One-Time Preprocessing

**IMPORTANT:** You must run preprocessing once to create the indexed GTF files:

```bash
bash preprocess_sqanti_tabix.sh
```

This script will:
1. Create tabix-indexed version of GENCODE GTF (2.9 GB -> 140 MB compressed)
2. Create tabix-indexed version of SQANTI GTF (1.0 GB -> 41 MB compressed)
3. Generate tabix indices for fast coordinate lookups

**Time:** ~10-15 minutes (one-time setup)

**Output files:**
- `/path/to/gencode.v49.primary_assembly.annotation.chrnamesedited.sorted.gtf.gz` (+ `.tbi` index)
- `/path/to/sqanti3_corrected.sorted.gtf.gz` (+ `.tbi` index)

### Required Files

Before running the annotation pipeline, ensure these files are accessible and paths are correctly set in `config.R` (see [Configuration](#configuration)):

| File Type | Description | Typical Size | Path in Script |
|-----------|-------------|--------------|----------------|
| **GENCODE GTF (indexed)** | Tabix-indexed GENCODE GTF | ~140 MB | `GENCODE_GTF_INDEXED` |
| **GENCODE GTF index** | Tabix index for GENCODE | ~500 KB | `<GENCODE_GTF_INDEXED>.tbi` |
| **GENCODE FASTA** | Reference transcript sequences | ~135 MB | `GENCODE_FASTA` |
| **SQANTI Classification** | PacBio isoform metadata | ~770 MB | `SQANTI_CLASSIFICATION` |
| **SQANTI GTF (indexed)** | Tabix-indexed PacBio isoforms | ~41 MB | `SQANTI_GTF_INDEXED` |
| **SQANTI GTF index** | Tabix index for SQANTI | ~224 KB | `<SQANTI_GTF_INDEXED>.tbi` |
| **SQANTI FASTA** | PacBio transcript sequences | ~497 MB | `SQANTI_FASTA` |
| **SQANTI Protein FASTA** | Pre-translated PacBio ORFs | varies | `SQANTI_PROTEIN_FASTA` |
| **DGEList RDS** | Expression data with sample metadata | ~29 MB | `DGELIST_RDS` |

---

## Isoform Annotation Pipeline

The pipeline extracts and integrates comprehensive transcript-level information for a gene of interest, combining:

1. **Reference annotations** from GENCODE (canonical transcripts)
2. **Novel isoforms** from PacBio long-read sequencing (SQANTI)
3. **Structural features** (exons, CDS regions, UTRs, splice junctions)
4. **Expression data** (configurable stratification across sample groups)
5. **Protein analysis** (CDS translation, UniProt comparison)

### Quick Start

```bash
# 1. Run for a gene of interest (TSV output only)
Rscript gene_isoform_annotation.R MTCL1

# 2. Include GTF and FASTA outputs
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta

# 3. Include protein translation and UniProt comparison
Rscript gene_isoform_annotation.R MTCL1 --protein

# 4. All outputs
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta --protein

# 5. Skip expression calculation (faster)
Rscript gene_isoform_annotation.R MTCL1 --no-expr

# 6. Use Ensembl ID instead of gene name
Rscript gene_isoform_annotation.R ENSG00000168502 --gtf --fasta --protein
```

**Output:** Tab-separated file with comprehensive isoform annotations (exon structures, CDS positions, UTR lengths, expression levels, protein analysis, and more).

**Performance:** ~30-60 seconds per gene using tabix-indexed GTF files.

**Prerequisite:** Run `bash preprocess_sqanti_tabix.sh` once to create indexed GTF files (~10-15 minutes).

### Command-Line Options

```bash
Rscript gene_isoform_annotation.R <GENE_NAME_OR_ID> [OPTIONS]
```

**Required Argument:**
- `GENE_NAME_OR_ID` - Gene symbol (e.g., `MTCL1`) or Ensembl ID (e.g., `ENSG00000168502`)

**Optional Flags:**
- `--gtf` - Generate GTF file with all genomic features
- `--fasta` - Generate FASTA file with transcript sequences
- `--protein` - Generate protein FASTA with CDS translations and UniProt comparison columns
- `--uniprot_file <path>` - Use a local UniProt FASTA file instead of fetching from the API
- `--no-expr` - Skip expression calculation (faster, no DGEList needed)

#### Examples

```bash
# Basic usage - TSV output only
Rscript gene_isoform_annotation.R MTCL1

# Include GTF for genome browser visualization
Rscript gene_isoform_annotation.R MTCL1 --gtf

# Include FASTA for sequence analysis
Rscript gene_isoform_annotation.R MTCL1 --fasta

# Protein translation and UniProt comparison
Rscript gene_isoform_annotation.R MTCL1 --protein

# Use a local UniProt FASTA instead of API fetch
Rscript gene_isoform_annotation.R MTCL1 --protein --uniprot_file /path/to/uniprot.fasta

# All outputs
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta --protein

# Structural annotations only (skip expression)
Rscript gene_isoform_annotation.R MTCL1 --no-expr

# Using Ensembl ID instead of gene name
Rscript gene_isoform_annotation.R ENSG00000168502

# Batch processing multiple genes
for gene in MTCL1 TP53 EGFR CFTR; do
  Rscript gene_isoform_annotation.R $gene --gtf --fasta --protein
done
```

#### Performance

**Runtime:** ~30-60 seconds per gene (without `--protein`); ~60-90 seconds with `--protein`

The pipeline uses tabix-indexed GTF files for fast coordinate-based queries:
- GENCODE extraction: <1 second (via tabix index)
- SQANTI extraction: <1 second (via tabix index)
- Expression calculation: ~15 seconds (DGEList loading and CPM calculation)
- Protein translation: ~15-30 seconds (FASTA loading + translation)
- Total: 30-90 seconds per gene

### Expression Configuration

Expression calculation is controlled by settings in `config.R`:

- `TRANSCRIPT_ID_COLUMN` - Column in `DGEList$genes` containing transcript IDs (default: `"txid"`)
- `STRATIFY_BY` - Columns in `DGEList$samples` to stratify expression by:
  - `NULL` → no stratification: `expr_mean`, `total_reads_mean`
  - `c("treatment")` → one-way: `expr_DMSO`, `expr_Smg1i`, ...
  - `c("treatment", "ct")` → two-way: `expr_DMSO_ddali`, `expr_Smg1i_at2`, ...
- `LEVEL_LABELS` - Optional renaming of factor levels for cleaner column names

**Process:**
1. Load DGEList and calculate normalized CPM (counts per million)
2. For each sample group (defined by `STRATIFY_BY` columns):
   - Calculate mean CPM across samples in the group
   - Calculate sum of raw read counts across samples in the group
3. Match isoforms by transcript ID (uses `TRANSCRIPT_ID_COLUMN`)
4. Round CPM values to 2 significant figures for readability

**Result:** Per-group expression and read count columns
- `expr_{group}` - Mean TMM-normalized CPM for each sample group
- `total_reads_{group}` - Sum of raw read counts for each sample group
- `total_reads_all_samples` - Grand total across all samples
- Missing values (NA) indicate isoform not detected in sequencing

### Output Files

#### TSV File (Primary Output)

**File:** `isoform_annotation_<GENE>.tsv`

Tab-separated table with one row per isoform. See `column_definitions.tsv` for the complete machine-readable column catalog.

##### Structural Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `isoform_id` | String | Transcript identifier | ENST00000306329.16 |
| `isoform_name` | String | Transcript name | MTCL1-201 |
| `source` | String | GENCODE or SQANTI | GENCODE |
| `chromosome` | String | Chromosome location | 18 |
| `strand` | String | Strand | + |
| `tss` | Integer | Transcription start site (strand-aware) | 8706502 |
| `tes` | Integer | Transcription end site (strand-aware) | 8832780 |
| `transcript_length` | Integer | Total length (bp, genomic span) | 126279 |
| `n_exons` | Integer | Number of exons | 14 |
| `exons` | String | Exon coords with CDS/UTR types | 8706502_8706713_+:CDS;... |
| `junctions` | String | Splice junction coords | 8706713_8718424_+;... |

##### CDS and Protein Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `cds_start` | Integer | CDS genomic start (always smaller coordinate) | 8706607 |
| `cds_end` | Integer | CDS genomic end (always larger coordinate) | 8831837 |
| `cds_length_nt` | Integer | CDS length in nucleotides | 5718 |
| `cds_length_aa` | Integer | CDS length in amino acids (floor of nt / 3) | 1906 |
| `in_frame` | Integer | CDS divisible by 3: 1 = yes, 0 = no, NA = no CDS | 1 |
| `utr5_length` | Integer | 5' UTR length (bp, strand-aware, genomic span) | 105 |
| `utr3_length` | Integer | 3' UTR length (bp, strand-aware, genomic span) | 943 |
| `is_protein_coding` | Boolean | Protein coding status | TRUE |
| `coding_source` | String | Source of coding annotation | GENCODE_biotype |

##### Expression Columns (per sample group)

| Column Pattern | Type | Description | Example |
|----------------|------|-------------|---------|
| `expr_{group}` | Float | Mean TMM-normalized CPM for sample group | 12.5 |
| `total_reads_{group}` | Integer | Sum of raw read counts for sample group | 347 |
| `total_reads_all_samples` | Integer | Grand total raw reads across all samples | 1523 |

Actual column names depend on `STRATIFY_BY` and `LEVEL_LABELS` in `config.R`. With the NMD defaults (`STRATIFY_BY = c("treatment", "ct")`), groups are `DMSO_ddali`, `Smg1i_ddali`, `DMSO_dd`, etc.

##### Protein Comparison Columns (requires `--protein`)

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `frame_matches_uniprot` | Boolean | Exact match or perfect subset of UniProt canonical | TRUE |
| `match_pct_uniprot` | Float | Percentage of 20-aa windows found in UniProt canonical | 93.4 |
| `matched_range_uniprot` | String | UniProt aa range of longest contiguous identical stretch | 781-1893 |

**Exon Format:** `start_end_strand:type`
- Types: `CDS`, `5UTR`, `3UTR`, `mixed`, `exon`
- Multiple exons separated by semicolons

**Junction Format:** `exon1_end_exon2_start_strand`
- Represents splice sites between consecutive exons

**Expression Values:**
- `NA` - Isoform not detected in sequencing
- `0` - Detected but not expressed
- `>0` - Expressed; CPM values rounded to 2 sig figs

#### GTF File (Optional)

**File:** `isoform_annotation_<GENE>.gtf`

Standard GTF format compatible with genome browsers:
- Gene, transcript, exon, CDS, and UTR features
- Both GENCODE and SQANTI isoforms (filtered to target gene only)
- Can be loaded into IGV, UCSC Genome Browser, etc.

#### FASTA File (Optional)

**File:** `isoform_annotation_<GENE>.fasta`

Standard FASTA format with transcript sequences:
```
>ENST00000306329.16
ATGGCGGCGGCG...
>PB.14156.3
ATGGAGGAGCTG...
```
- Contains all isoforms (GENCODE + SQANTI)
- Headers are simple transcript IDs
- Extracted from two sources (GENCODE and SQANTI FASTA files)

#### Protein FASTA File (Optional)

**File:** `protein_sequences_<GENE>.fasta`

Amino acid sequences for all translated isoforms:
```
>sp|Q9Y4B5|MTCL1_HUMAN Microtubule cross-linking factor 1 OS=Homo sapiens ...
MEEAVLQAPGP...
>ENST00000306329.16
MEEAVLQAPGP...
>PB.14156.3
MFHALLQDSGR...
```
- First entry is the UniProt canonical (with full Swiss-Prot header)
- Subsequent entries are CDS translations keyed by isoform ID
- Translated proteins may include a trailing `*` (stop codon)

### Sample Output

#### Example: MTCL1 Gene

**Input:**
```bash
Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta --protein
```

**Output Files:**
```
isoform_annotation_MTCL1.tsv       # 58 isoforms with full annotation columns
isoform_annotation_MTCL1.gtf       # 640 genomic features
isoform_annotation_MTCL1.fasta     # 58 transcript sequences
protein_sequences_MTCL1.fasta      # 42 protein sequences (1 UniProt + 41 translated)
```

**Summary Statistics:**
- 21 GENCODE reference isoforms
- 37 SQANTI novel isoforms
- 58 total isoforms annotated
- 41 protein coding, 17 non-coding
- 8-9 isoforms with detectable expression per cell type
- 3 isoforms with exact/subset match to UniProt canonical protein

---

## BAM Visualization & IGV Sessions

Tools for visualizing genome-aligned (hg38) PacBio long-read BAMs in IGV and as sashimi plots. The `extract_gene_region.sh` script extracts gene regions from genome-aligned BAMs and produces small BAM slices with an IGV session file.

### Quick Start

```bash
# Extract a gene region from genome-aligned BAMs
bash extract_gene_region.sh --region chr18:8705271-8832780 \
    --bam-dir ../data/bams/genome_hg38/ \
    --output-dir ./igv_genome/ \
    --gtf isoform_annotation_MTCL1.gtf

# Open the resulting session in IGV
# File > Open Session > igv_genome/igv_session.xml
```

**Dependencies:** `samtools` (required), `tabix` (only for `--gene` lookup)

### extract_gene_region.sh

A portable shell script that extracts a gene region from genome-aligned BAMs into small, self-contained files suitable for IGV visualization.

**Usage:**

```bash
# Extract by genomic coordinates (primary interface)
bash extract_gene_region.sh --region chr18:8705271-8832780 \
    --bam-dir ../data/bams/genome_hg38/ \
    --output-dir ./igv_genome/ \
    --gtf isoform_annotation_MTCL1.gtf

# Extract specific BAM files
bash extract_gene_region.sh --region chr18:8705271-8832780 \
    --bam /path/to/sample1.bam --bam /path/to/sample2.bam \
    --output-dir ./igv_genome/

# Look up gene coordinates from GENCODE GTF (slower, greps full file)
bash extract_gene_region.sh --gene MTCL1 \
    --gencode-gtf /path/to/gencode.v49.annotation.gtf.gz \
    --bam-dir ../data/bams/genome_hg38/ \
    --output-dir ./igv_genome/

# Adjust padding (default: 5000 bp on each side)
bash extract_gene_region.sh --region chr18:8705271-8832780 \
    --pad 10000 --bam-dir ../data/bams/genome_hg38/ \
    --output-dir ./igv_genome/
```

**Options:**

| Flag | Description |
|------|-------------|
| `--region CHR:START-END` | Genomic region (e.g., `chr18:8705271-8832780`) |
| `--gene GENE_NAME` | Gene name lookup (requires `--gencode-gtf`) |
| `--bam-dir DIR` | Directory containing `.bam` files |
| `--bam FILE` | Specific BAM file (repeatable) |
| `--output-dir DIR` | Output directory (required) |
| `--gtf FILE` | Gene annotation GTF to include in IGV session |
| `--gencode-gtf FILE` | Bgzipped GENCODE GTF for `--gene` lookup |
| `--pad N` | Padding around region in bp (default: 5000) |
| `--display-mode MODE` | IGV display mode: SQUISHED, EXPANDED, COLLAPSED (default: SQUISHED) |

**Features:**
- Detects chromosome naming convention (`chr18` vs `18`) in BAM headers and adjusts queries automatically
- Validates BAM indexes; creates them if missing
- Warns when a BAM has 0 reads in the extracted region
- Generates IGV session XML with tracks sorted by cell type, donor, treatment (DMSO before Smg1i)
- Produces small BAM slices (~30-50 KB per sample for a single gene, down from ~3-4 GB)

**Output:**

```
igv_genome/
├── Sample13_DD_017Q_DMSO.aligned.bam       # Extracted BAM slice
├── Sample13_DD_017Q_DMSO.aligned.bam.bai   # BAM index
├── Sample14_DD_017Q_Smg1i.aligned.bam
├── Sample14_DD_017Q_Smg1i.aligned.bam.bai
├── gene_annotation.gtf                      # Copied from --gtf
└── igv_session.xml                          # Open in IGV
```

### Sashimi Plots with ggsashimi

For publication-quality splice junction visualization with intron shrinking, use [ggsashimi](https://github.com/guigolab/ggsashimi) via Docker:

```bash
# Pull the Docker image (one-time)
docker pull guigolab/ggsashimi

# Run with intron shrinking, colored by treatment
docker run --platform linux/amd64 \
  -v "$(pwd)/igv_genome:/data" \
  -w /data \
  guigolab/ggsashimi \
  -b /data/input_bams.tsv \
  -c chr18:8700271-8837780 \
  -g /data/gene_annotation_filtered.gtf \
  --shrink \
  -C 3 -M 2 --fix-y-scale \
  --height 4 --ann-height 6 --width 20 --base-size 16 --alpha 0.9 \
  -P /data/palette.tsv \
  -o /data/sashimi_MTCL1 \
  -F pdf
```

**Input files** (in `igv_genome/`):

- `input_bams.tsv` — Tab-separated: `label<TAB>/data/sample.bam<TAB>group` (group is used for coloring with `-C 3`)
- `palette.tsv` — Tab-separated: `hex_color<TAB>group_name` (e.g., `#2166ac<TAB>DMSO`)
- `gene_annotation_filtered.gtf` — GTF filtered to expressed isoforms only (reduces annotation clutter)

**Aggressive intron shrinking:** ggsashimi's default `--shrink` uses exponent 0.7. For genes with large introns, a patched version (`ggsashimi_patched.py`) with exponent 0.3 produces more compact plots. Run via `--entrypoint python` with the patched script mounted.

**Filtering the annotation GTF:** The full annotation GTF contains 64 isoforms (43 PacBio + 21 GENCODE), which overwhelms the annotation panel. Filter to expressed isoforms using primary alignment read counts (see `igv_genome/expressed_isoforms.txt`).

---

## Reference

### How the Annotation Pipeline Works

#### Step 1: Gene Identification
- Accepts either gene name (e.g., "MTCL1") or Ensembl ID (e.g., "ENSG00000168502")
- Searches GENCODE GTF to find gene ID and confirm canonical gene name
- Extracts gene coordinates (chromosome, start, end) for tabix queries
- Validates gene exists in reference genome

#### Step 2: GENCODE Isoform Extraction
- Queries tabix-indexed GENCODE GTF by gene coordinates (<1 second)
- Filters results to the target gene ID (excludes overlapping genes)
- Parses GTF to extract transcript features
- For each transcript:
  - **Exon coordinates** with functional types (CDS, 5'UTR, 3'UTR, mixed)
  - **Splice junctions** between consecutive exons
  - **TSS/TES** (transcription start/end sites, strand-aware)
  - **CDS boundaries** (genomic start/end of coding sequence)
  - **UTR lengths** calculated from TSS->CDS (5'UTR) and CDS->TES (3'UTR)
  - **Protein coding status** from transcript biotype and CDS presence
  - **CDS length** in nucleotides and amino acids
  - **In-frame status** (whether CDS length is divisible by 3)

**Key feature:** Strand-aware calculations ensure correct UTR assignment:
- **Plus strand** (+): TSS < CDS_start < CDS_end < TES
  - 5'UTR = CDS_start - TSS
  - 3'UTR = TES - CDS_end
- **Minus strand** (-): TES < CDS_start < CDS_end < TSS
  - 5'UTR = TSS - CDS_end
  - 3'UTR = CDS_start - TES

#### Step 3: SQANTI Isoform Extraction

**Tabix-indexed coordinate queries (<1 second):**
- Queries preprocessed, coordinate-indexed SQANTI GTF
- Retrieves only features overlapping gene coordinates
- Uses same coordinates as GENCODE for consistent coverage

For each SQANTI isoform:
- Extracts exon structures (simpler than GENCODE, no UTR subtypes)
- Gets coding status from SQANTI classification
- Calculates CDS positions and UTRs if available
- Identifies novel isoforms (PB.* IDs) not in reference

#### Step 4: Expression Calculation (Optional)

**Data Source:** DGEList object containing:
- TMM-normalized counts from edgeR
- Sample metadata with stratification columns (configurable in `config.R`)

See [Expression Configuration](#expression-configuration) for details on `STRATIFY_BY` and related settings.

#### Step 5: Protein Translation and UniProt Comparison (Optional, `--protein`)

When the `--protein` flag is used, the pipeline performs CDS translation and compares each predicted protein to the UniProt canonical sequence. See [Protein Analysis Algorithms](#protein-analysis-algorithms) for full details.

**Process:**
1. Fetch the canonical UniProt protein sequence (via REST API or local file)
2. Load transcript FASTA sequences (GENCODE + SQANTI)
3. Translate each protein-coding isoform's CDS to amino acids
4. Compare each translated protein to the UniProt canonical
5. Write all protein sequences to a FASTA file

**Translation sources (in order of priority):**
- **GENCODE isoforms**: CDS features from GTF mapped to transcript coordinates, then translated via `Biostrings::translate()`
- **SQANTI isoforms**: CDS features from SQANTI GTF if available, otherwise falls back to pre-translated ORFs from `sqanti3_corrected.faa`

**UniProt comparison columns added to TSV:**
- `frame_matches_uniprot` - TRUE if exact match or perfect substring of canonical
- `match_pct_uniprot` - Percentage of sequence matching the canonical (sliding window)
- `matched_range_uniprot` - UniProt amino acid range of the longest identical stretch

#### Step 6: Output Generation

**TSV (always created):**
- One row per isoform with all annotation columns
- Expression values rounded to 2 significant figures
- Protein comparison columns included when `--protein` is used
- Ready for downstream analysis in R/Python

**GTF (optional, `--gtf` flag):**
- Contains all genomic features for visualization
- Combines GENCODE + SQANTI features (filtered to target gene only)
- Compatible with genome browsers (IGV, UCSC)

**FASTA (optional, `--fasta` flag):**
- Transcript sequences for all isoforms
- Extracts from two sources:
  - GENCODE transcripts -> GENCODE FASTA
  - SQANTI transcripts -> SQANTI FASTA
- Handles pipe-delimited GENCODE headers automatically
- Useful for sequence analysis, primer design, structure prediction

**Protein FASTA (optional, `--protein` flag):**
- Protein sequences for all translated isoforms
- First entry is the UniProt canonical sequence (with full header)
- Subsequent entries are CDS translations keyed by isoform ID
- Useful for multiple sequence alignment, domain analysis, structure prediction

### Protein Analysis Algorithms

#### CDS-to-Transcript Coordinate Mapping

To translate a CDS, the pipeline must map genomic CDS coordinates to positions within the mRNA transcript sequence. This is non-trivial because CDS regions may span multiple exons separated by introns.

**Algorithm (`get_cds_transcript_coords`):**

1. Sort exons by genomic position
2. Determine transcript order: same as genomic for plus strand, reversed for minus strand
3. Build cumulative exon length array in transcript order
4. For each exon (in transcript order), check if the CDS start or end boundary falls within it:
   - **Plus strand**: CDS start (lowest genomic coordinate) maps to transcript CDS start; CDS end (highest) maps to transcript CDS end
   - **Minus strand**: CDS end (highest genomic coordinate) maps to transcript CDS start (first codon); CDS start (lowest) maps to transcript CDS end (last codon)
5. Calculate the 1-based transcript-relative position: `cumulative_length_before_exon + offset_within_exon + 1`
6. **Sanity check**: Verify that the mapped CDS length (`tx_end - tx_start + 1`) equals the sum of all CDS feature widths from the GTF. If mismatched, the translation is skipped with a warning.

**Translation:**

The mapped CDS subsequence is extracted from the full transcript sequence using `Biostrings::subseq()` and translated with `Biostrings::translate()`. Isoforms with CDS lengths not divisible by 3 will fail translation (caught by `tryCatch`) and are flagged via the `in_frame` column.

#### UniProt Comparison: `frame_matches_uniprot`

For each translated protein, the pipeline checks whether it is an **exact match** or a **perfect substring** of the UniProt canonical sequence:

```
is_match = (isoform_protein == uniprot_canonical) OR
           (isoform_protein is a contiguous substring of uniprot_canonical)
```

The trailing stop codon (`*`) produced by `Biostrings::translate()` is stripped before comparison, since UniProt sequences do not include it.

**Values:**
- `TRUE` - The translated protein is identical to, or perfectly contained within, the UniProt canonical
- `FALSE` - The translated protein diverges from the canonical (frame shifts, alternative exons, etc.)
- `NA` - Not protein-coding, or translation failed

#### UniProt Comparison: `match_pct_uniprot`

Measures what fraction of the translated protein sequence can be found in the UniProt canonical, using a **sliding window approach**:

1. Slide a 20-amino-acid window across the translated protein, one position at a time
2. For each window position, check if that exact 20-aa substring exists anywhere in the UniProt canonical sequence (using exact string matching)
3. Calculate: `match_pct = 100 * (number of matching windows) / (total windows)`

Where `total windows = protein_length - 20 + 1`.

For proteins shorter than 20 aa, a simple full-sequence substring check is used instead.

**Interpretation:**
- `100%` - Every 20-aa window is found in the canonical (exact match or perfect subset)
- `90-99%` - Most of the protein matches; divergence at one or both termini (common with alternative first/last exons)
- `<80%` - Substantial sequence divergence, likely an alternative reading frame or heavily rearranged

#### UniProt Comparison: `matched_range_uniprot`

Reports the **UniProt amino acid coordinates** of the longest contiguous stretch of identical sequence between the isoform and the canonical:

1. Using the same 20-aa sliding window results, find the longest run of consecutive matching windows (via run-length encoding)
2. The contiguous matching stretch spans from the first position of the run to the last position of the run + 19 (window size - 1)
3. Look up where this stretch falls in the UniProt canonical sequence
4. Report as `start-end` (1-based UniProt coordinates)

**Example:** `781-1893` means the longest identical stretch maps to UniProt amino acids 781 through 1893 (1113 aa of continuous identity).

### Data Flow Diagram

```
USER INPUT
    |
    v
+-------------------------------------------+
|  Gene Name (MTCL1)                        |
|  OR Ensembl ID (ENSG00000168502)          |
|  + Optional Flags                         |
|    --gtf --fasta --protein --no-expr      |
+-------------------------------------------+
    |
    v
+---------------------------------------------------------------+
| STEP 1: GENE IDENTIFICATION                                    |
|   - Search GENCODE GTF (grep gene feature)                    |
|   - Resolve gene name <-> Ensembl ID                           |
|   - Extract gene coordinates: chr18:8705271-8832780           |
|   - Result: ENSG00000168502.19 -> "MTCL1"                      |
+---------------------------------------------------------------+
    |
    v
+---------------------------------------------------------------+
| STEP 2: GENCODE EXTRACTION                                     |
|   Input: Tabix-indexed GENCODE GTF (140 MB compressed)        |
|   Method: tabix query by coordinates (<1 second)              |
|   - Parse exons (with CDS/UTR types)                          |
|   - Calculate junctions, TSS/TES, CDS positions               |
|   - Calculate UTR lengths (strand-aware)                       |
|   - Determine protein coding status and in-frame status        |
|   Output: 21 GENCODE isoforms with full annotations           |
+---------------------------------------------------------------+
    |
    v
+---------------------------------------------------------------+
| STEP 3: SQANTI EXTRACTION                                      |
|   Input: Tabix-indexed SQANTI GTF + Classification            |
|   Method: tabix query by coordinates (<1 second)              |
|   - Extract exon positions from GTF                           |
|   - Get coding status from SQANTI classification              |
|   - Calculate CDS positions and UTRs                          |
|   Output: 37 SQANTI isoforms with annotations                 |
+---------------------------------------------------------------+
    |
    v
+---------------------------------------------------------------+
| STEP 4: EXPRESSION CALCULATION (if --no-expr NOT set)         |
|   Input: DGEList RDS + STRATIFY_BY config                      |
|   - Load TMM-normalized CPM + raw counts                       |
|   - For each sample group (defined by STRATIFY_BY):            |
|     - Mean CPM across samples in the group                     |
|     - Sum of raw read counts across samples in the group       |
|   - Round CPM to 2 significant figures                         |
|   Output: expr/total_reads per group + total_reads_all_samples |
+---------------------------------------------------------------+
    |
    v
+---------------------------------------------------------------+
| STEP 5: PROTEIN ANALYSIS (if --protein set)                   |
|   - Fetch UniProt canonical sequence (API or local file)       |
|   - Translate CDS for all protein-coding isoforms:             |
|     - GENCODE: GTF CDS -> transcript coords -> translate       |
|     - SQANTI: GTF CDS or sqanti3_corrected.faa fallback        |
|   - Compare each protein to UniProt canonical:                 |
|     - Exact/subset match (frame_matches_uniprot)               |
|     - Sliding 20-aa window match % (match_pct_uniprot)         |
|     - Longest contiguous match range (matched_range_uniprot)   |
|   Output: protein FASTA + 3 comparison columns in TSV          |
+---------------------------------------------------------------+
    |
    v
+---------------------------------------------------------------+
| STEP 6: OUTPUT GENERATION                                      |
|                                                                |
|  TSV (always)   GTF (--gtf)   FASTA (--fasta)  Protein FASTA  |
|  58 isoforms    640 features  58 sequences     42 proteins     |
|  All columns    For genome    For sequence     1 UniProt +     |
|  Expression +   browsers      analysis         41 translated   |
|  protein data                                                  |
+---------------------------------------------------------------+
    |
    v
FINAL OUTPUT: isoform_annotation_<GENE>.tsv
              + isoform_annotation_<GENE>.gtf
              + isoform_annotation_<GENE>.fasta
              + protein_sequences_<GENE>.fasta
```

### NMD Study Terminology

**DMSO**: Dimethyl sulfoxide, used as vehicle control in NMD inhibitor treatment. DMSO expression represents **baseline expression levels** in untreated cells (before NMD pathway inhibition).

**NMD**: Nonsense-Mediated Decay, a cellular quality control pathway that degrades aberrant transcripts. Inhibiting NMD (with Smg1i) reveals transcripts normally targeted for degradation.

**Isoform**: Alternative transcript variant of a gene, created through alternative splicing, alternative TSS/TES, or other mechanisms.

#### Cell Types (Primary Human Lung Cells)

| Code | Full Name | Description |
|------|-----------|-------------|
| **DD** | Basal large airway epithelial cells | Differentiated cells from large airways |
| **DD_ALI** | Basal airway epithelial cells (ALI) | Basal airway cells cultured at air-liquid interface |
| **DO_ALI** | Basal small airway cells (ALI) | Small airway cells cultured at air-liquid interface |
| **AT2** | Alveolar epithelial type 2 cells | Surfactant-producing cells in alveoli |
| **FB** | Lung fibroblasts | Connective tissue cells |
| **MV** | Lung microvascular cells | Endothelial cells from lung microvasculature |

**ALI (Air-Liquid Interface)**: Culture condition where cells are grown at the interface between air and liquid medium, promoting differentiation and mucus production.

#### Data Sources

**GENCODE**: Reference human genome annotation consortium providing high-quality gene and transcript annotations.

**SQANTI**: Software for characterizing PacBio long-read transcripts, including novel isoforms not in reference annotations.

**UniProt**: The Universal Protein Resource, providing reviewed (Swiss-Prot) canonical protein sequences.

---

## Troubleshooting

### Common Issues

**Error: "Gene not found in GENCODE"**
```
Solution:
- Check gene name spelling (case-sensitive)
- Try Ensembl ID: Rscript gene_isoform_annotation.R ENSG00000168502
- Verify gene is in GENCODE v49
```

**Error: "DGEList RDS not found" or column not found in DGEList**
```
Solution:
- If you don't have expression data, use --no-expr flag
- Rscript gene_isoform_annotation.R MTCL1 --no-expr
- Check TRANSCRIPT_ID_COLUMN and STRATIFY_BY in config.R match your DGEList
```

**Error: "Missing required files" or "tabix index not found"**
```
Solution:
- Run preprocessing first: bash preprocess_sqanti_tabix.sh
- This creates the required tabix-indexed GTF files
- Takes 10-15 minutes but only needed once
```

**Error: "tabix command not found" or "bgzip command not found"**
```
Solution:
- Install htslib: conda install -c bioconda htslib
```

**Warning: "N isoforms not found in DGEList"**
```
This is normal:
- Newer GENCODE versions may not be in DGE analysis
- Novel SQANTI isoforms may not be quantified
- Affected isoforms will have NA expression values
```

**Error: Package 'Biostrings' not found (with --fasta or --protein)**
```
Solution:
- Install Biostrings: BiocManager::install("Biostrings")
- OR run without --fasta / --protein flags
```

**Warning: "last base was ignored" during protein translation**
```
This is normal:
- One or more isoforms have a CDS length not divisible by 3
- The in_frame column flags these (in_frame = 0)
- Translation proceeds using complete codons only
```

### Performance Issues

**Script taking longer than expected?**
1. Check if tabix preprocessing was done: `ls sqanti3_corrected.sorted.gtf.gz.tbi`
2. Monitor system resources: `top` or `htop`
3. Large genes (>50 isoforms) take longer
4. Expression calculation adds ~15-30 seconds - use `--no-expr` if not needed
5. Protein translation adds ~15-30 seconds (FASTA file loading)

**Out of memory errors?**
1. Close other applications
2. Use `--no-expr` to reduce memory usage
3. Process genes one at a time instead of batch mode

---

## Citation

If you use this pipeline in your research, please cite:

- **GENCODE**: Frankish et al. (2021). GENCODE 2021. *Nucleic Acids Research*, 49(D1), D916-D923.
- **SQANTI3**: Tardaguila et al. (2018). SQANTI: extensive characterization of long-read transcript sequences. *Genome Research*, 28(3), 396-411.
- **edgeR**: Robinson et al. (2010). edgeR: a Bioconductor package for differential expression analysis. *Bioinformatics*, 26(1), 139-140.
- **UniProt**: The UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research*, 51(D1), D523-D531.

---

## Version History

- **v3.2** (2026-02-27)
  - Added `environment.yml` for reproducible conda setup (R, all R packages, htslib, samtools)
  - Removed `check_compiler()` from R script (irrelevant with conda binaries)
  - README now offers two setup paths: conda environment (recommended) vs manual installation

- **v3.1** (2026-02-27)
  - Reorganized README for usability: added Table of Contents, promoted BAM visualization to peer section, consolidated deep-dive material into Reference section, moved Prerequisites/Installation up into Setup

- **v3.0** (2026-02-26)
  - Standalone repository (`lr_igv_from_bam`) — no longer embedded in parent project
  - Config-based paths: `config.R` and `config.sh` replace all hardcoded paths
  - Template files (`config.example.R`, `config.example.sh`) for easy setup
  - Removed legacy `gene_isoform_annotation_fast.R`

- **v2.1** (2026-02-26)
  - Added `extract_gene_region.sh` for genome-aligned BAM extraction and IGV session generation
  - Added ggsashimi sashimi plot workflow (Docker-based, with intron shrinking)
  - Expression-filtered annotation GTF for cleaner sashimi plot annotation panels

- **v2.0** (2026-02-20)
  - Added `--protein` flag for CDS translation and protein FASTA output
  - Added `--uniprot_file` option for local UniProt FASTA input
  - Added UniProt comparison columns: `frame_matches_uniprot`, `match_pct_uniprot`, `matched_range_uniprot`
  - Renamed `cds_length` to `cds_length_nt` and `protein_length` to `cds_length_aa`
  - Added `in_frame` column (CDS divisible by 3)
  - Added Smg1i expression and raw read count columns alongside DMSO
  - Added `total_reads_all_samples` summary column
  - Added `column_definitions.tsv` machine-readable column catalog
  - Improved SQANTI protein FASTA header parsing (tab-delimited)
  - Performance: list accumulation replaces O(n^2) rbind loops
  - Robustness: CDS length sanity check, URL encoding, API timeout, gene-filtered GTF output

- **v1.0** (2026-02-14)
  - Initial production release
  - Command-line interface with gene name or Ensembl ID
  - Optional GTF and FASTA outputs via flags
  - Tabix optimization with automatic fallback
  - DMSO-specific expression from DGEList
  - CDS positions and strand-aware UTR calculations
  - Expression values rounded to 2 significant figures
  - Comprehensive error checking and diagnostics

---

## Support

For questions or issues:
1. Check the [Troubleshooting](#troubleshooting) section above
2. Verify file paths in `config.R` (see [Configuration](#configuration))
3. Ensure all required files are accessible
4. Review diagnostic messages in output

---

## License

This pipeline was developed as part of the NMD research project. Please cite appropriately if used in publications.
