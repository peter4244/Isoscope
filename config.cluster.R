# Gene Isoform Annotation Pipeline - Configuration
# Cluster paths — rename to config.R:
#   mv config.cluster.R config.R

# GENCODE reference files (gene lookup file ships with the repo — no config needed)
GENCODE_GTF_INDEXED   <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GENCODE_FASTA         <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/reference_files/gencode.v49.transcripts.fa.gz"

# SQANTI PacBio files
SQANTI_CLASSIFICATION <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/code/sqanti3_by_chromosome/projects/nmd_lungcells/results/nmd_lungcells_classification.txt"
SQANTI_GTF_INDEXED    <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/code/sqanti3_by_chromosome/projects/nmd_lungcells/results/nmd_lungcells_corrected.sorted.gtf.gz"
SQANTI_FASTA          <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/code/sqanti3_by_chromosome/projects/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
SQANTI_PROTEIN_FASTA  <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/code/sqanti3_by_chromosome/projects/nmd_lungcells/results/nmd_lungcells_corrected.faa"

# Expression data
DGELIST_RDS           <- "/proj/regeps/regep00/studies/ExternalCellLines/analyses/repjc/Randell_Lung_Cells_2025/results/isocall_dge/dgelist_isocall_2026.3.1.rds"

# === Expression settings (used by Step 4, ignored with --no-expr) ===

# Column in DGEList$genes containing transcript IDs
TRANSCRIPT_ID_COLUMN <- "transcript_id"

# Columns in DGEList$samples to stratify expression by.
#   NULL            -> no stratification: expr_mean, total_reads
#   c("treatment")  -> one-way: expr_DMSO, expr_Smg1i, ...
#   c("treatment", "ct") -> two-way: expr_DMSO_ddali, expr_Smg1i_at2, ...
STRATIFY_BY <- c("treatment", "ct")

# Optional: rename factor levels for cleaner column names.
# Named list of named character vectors, keyed by column name.
# NULL for a column = use raw levels as-is.
# If provided for a column, ALL levels in that column must have a mapping.
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
