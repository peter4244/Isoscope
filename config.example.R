# Gene Isoform Annotation Pipeline - Configuration
# Copy this file to config.R and edit paths for your system:
#   cp config.example.R config.R
#
# All paths must be absolute. See README.md for file descriptions.

# GENCODE reference files (gene lookup file ships with the repo — no config needed)
GENCODE_GTF_INDEXED   <- "/path/to/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GENCODE_FASTA         <- "/path/to/gencode.v49.transcripts.fa.gz"

# SQANTI PacBio files
SQANTI_CLASSIFICATION <- "/path/to/nmd_lungcells_classification.txt"
SQANTI_GTF_INDEXED    <- "/path/to/nmd_lungcells_corrected.sorted.gtf.gz"
SQANTI_FASTA          <- "/path/to/nmd_lungcells_corrected.fasta.gz"
SQANTI_PROTEIN_FASTA  <- "/path/to/nmd_lungcells_corrected.faa"

# Expression data (not needed with --no-expr)
DGELIST_RDS           <- "/path/to/dge_isoform.rds"

# === Expression settings (used by Step 4, ignored with --no-expr) ===

# Column in DGEList$genes containing transcript IDs
TRANSCRIPT_ID_COLUMN <- "transcript_id"

# Columns in DGEList$samples to stratify expression by.
#   NULL            -> no stratification: expr_mean, total_reads
#   c("treatment")  -> one-way: expr_DMSO, expr_Smg1i, ...
#   c("treatment", "ct") -> two-way: expr_DMSO_ddali, expr_Smg1i_at2, ...
STRATIFY_BY <- c("treatment", "ct")

# Optional: rename factor levels for cleaner column names and control order.
# Named list of named character vectors, keyed by column name.
# NULL for a column = use factor levels() or unique() order as-is.
# If provided for a column, ALL levels in that column must have a mapping.
# The order of entries controls the output column order.
LEVEL_LABELS <- list(
  treatment = c("DMSO" = "DMSO", "Smg1i" = "Smg1i"),
  ct = c(
    "DD_ALI" = "ddali",
    "DD"     = "dd",
    "DO"     = "doali",
    "AT"     = "at2",
    "FB"     = "fb",
    "MV"     = "mv"
  )
)
