# Gene Isoform Annotation Pipeline - Configuration
# Copy this file to config.R and edit paths for your system:
#   cp config.example.R config.R
#
# All paths must be absolute. See README.md for file descriptions.

# GENCODE reference files
GENCODE_GTF_UNINDEXED <- "/path/to/gencode.v49.primary_assembly.annotation.chrnamesedited.gtf"
GENCODE_GTF_INDEXED   <- "/path/to/gencode.v49.primary_assembly.annotation.chrnamesedited.sorted.gtf.gz"
GENCODE_FASTA         <- "/path/to/gencode.v49.transcripts.fa.gz"

# SQANTI PacBio files
SQANTI_CLASSIFICATION <- "/path/to/sqanti3_classification.txt"
SQANTI_GTF_INDEXED    <- "/path/to/sqanti3_corrected.sorted.gtf.gz"
SQANTI_FASTA          <- "/path/to/merged_collapsed.fa.gz"
SQANTI_PROTEIN_FASTA  <- "/path/to/sqanti3_corrected.faa"

# Expression data (not needed with --no-expr)
DGELIST_RDS           <- "/path/to/dge_isoform.rds"
