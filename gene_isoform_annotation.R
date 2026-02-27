#!/usr/bin/env Rscript
# Gene Isoform Annotation Pipeline - Production Version
# Extracts comprehensive isoform annotations for genes of interest
#
# Features:
# - Combines GENCODE reference annotations with SQANTI PacBio isoforms
# - Detailed exon structures with CDS/UTR annotations
# - Configurable expression stratification (from DGEList)
# - Fast processing using grep pre-filtering and tabix indexing
# - Optional GTF and FASTA outputs
#
# Usage:
#   Rscript gene_isoform_annotation.R <GENE_NAME_OR_ID>
#   Example: Rscript gene_isoform_annotation.R MTCL1
#   Example: Rscript gene_isoform_annotation.R ENSG00000168502
#
# Prerequisites:
#   1. Run preprocess_sqanti_tabix.sh once to create tabix-indexed SQANTI GTF
#   2. Ensure all reference files and DGEList RDS are available

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(edgeR)
})

# =============================================================================
# COMMAND-LINE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage: Rscript gene_isoform_annotation.R <GENE_NAME_OR_ID> [OPTIONS]\n")
  cat("\nArguments:\n")
  cat("  GENE_NAME_OR_ID    Gene name (e.g., MTCL1) or Ensembl ID (e.g., ENSG00000168502)\n")
  cat("\nOptions:\n")
  cat("  --gtf              Output GTF file with all isoforms\n")
  cat("  --fasta            Output FASTA file with isoform sequences\n")
  cat("  --protein          Output protein FASTA (translated CDS + Uniprot canonical)\n")
  cat("  --uniprot_file <f> Local Uniprot FASTA (default: fetch from API)\n")
  cat("  --no-expr          Skip expression calculation\n")
  cat("\nExamples:\n")
  cat("  Rscript gene_isoform_annotation.R MTCL1\n")
  cat("  Rscript gene_isoform_annotation.R MTCL1 --gtf --fasta --protein\n")
  cat("  Rscript gene_isoform_annotation.R ENSG00000168502 --no-expr\n")
  stop("No gene name or ID provided")
}

GENE_INPUT <- args[1]

# Parse optional flags
OUTPUT_GTF <- "--gtf" %in% args
OUTPUT_FASTA <- "--fasta" %in% args
OUTPUT_PROTEIN <- "--protein" %in% args
INCLUDE_EXPRESSION <- !"--no-expr" %in% args

# Parse --uniprot_file option (takes a value argument)
UNIPROT_FILE <- NULL
if ("--uniprot_file" %in% args) {
  uf_idx <- which(args == "--uniprot_file")
  if (uf_idx < length(args)) {
    UNIPROT_FILE <- args[uf_idx + 1]
  } else {
    stop("ERROR: --uniprot_file requires a file path argument")
  }
}

# Determine if input is Ensembl ID or gene name
if (grepl("^ENSG[0-9]+", GENE_INPUT)) {
  GENE_ID <- sub("\\..*", "", GENE_INPUT)  # Remove version if present
  GENE_NAME <- NULL  # Will be determined from GTF
  SEARCH_BY <- "id"
} else {
  GENE_NAME <- GENE_INPUT
  GENE_ID <- NULL  # Will be determined from GTF
  SEARCH_BY <- "name"
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# Note: OUTPUT_GTF, OUTPUT_FASTA, and INCLUDE_EXPRESSION are now set via command-line flags

# Load file paths from config (copy config.example.R -> config.R and edit)
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
config_path <- file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  stop("config.R not found. Copy config.example.R to config.R and edit paths for your system.")
}
source(config_path)

# Gene lookup file ships with the repo
GENCODE_GENE_LOOKUP <- file.path(script_dir, "gencode_v49_genes.gtf")

# Output files (will be set after gene name is resolved)
OUTPUT_TSV <- NULL
OUTPUT_GTF_FILE <- NULL
OUTPUT_FASTA_FILE <- NULL
OUTPUT_PROTEIN_FILE <- NULL

# =============================================================================
# HELPER FUNCTIONS FOR GENCODE PARSING
# =============================================================================

get_exon_structure_annotated <- function(gtf, transcript_id) {
  # Extract exon positions with CDS/UTR annotations
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]
  if (length(tx_exons) == 0) return(list(exons = NA, junctions = NA, n_exons = 0))

  tx_exons <- tx_exons[order(start(tx_exons))]
  strand_char <- as.character(strand(tx_exons)[1])

  # Get CDS/UTR annotations
  cds <- gtf[gtf$type == "CDS" & gtf$transcript_id == transcript_id]
  utr5 <- gtf[gtf$type == "five_prime_utr" & gtf$transcript_id == transcript_id]
  utr3 <- gtf[gtf$type == "three_prime_utr" & gtf$transcript_id == transcript_id]

  # Annotate each exon
  exon_types <- rep("exon", length(tx_exons))
  for (i in seq_along(tx_exons)) {
    exon <- tx_exons[i]
    cds_overlap <- any(overlapsAny(exon, cds))
    utr5_overlap <- any(overlapsAny(exon, utr5))
    utr3_overlap <- any(overlapsAny(exon, utr3))

    if (cds_overlap && (utr5_overlap || utr3_overlap)) {
      exon_types[i] <- "mixed"
    } else if (cds_overlap) {
      exon_types[i] <- "CDS"
    } else if (utr5_overlap) {
      exon_types[i] <- "5UTR"
    } else if (utr3_overlap) {
      exon_types[i] <- "3UTR"
    }
  }

  exon_names <- paste0(start(tx_exons), "_", end(tx_exons), "_", strand_char, ":", exon_types)

  # Junctions
  if (length(tx_exons) > 1) {
    junctions <- paste0(end(tx_exons[-length(tx_exons)]), "_", start(tx_exons[-1]), "_", strand_char)
    junctions_str <- paste(junctions, collapse = ";")
  } else {
    junctions_str <- NA
  }

  return(list(
    exons = paste(exon_names, collapse = ";"),
    junctions = junctions_str,
    n_exons = length(tx_exons)
  ))
}

get_exon_structure_simple <- function(gtf, transcript_id) {
  # Extract exon positions without CDS/UTR annotations (for SQANTI)
  tx_exons <- gtf[gtf$type == "exon" & gtf$transcript_id == transcript_id]
  if (length(tx_exons) == 0) return(list(exons = NA, junctions = NA, n_exons = 0))

  tx_exons <- tx_exons[order(start(tx_exons))]
  strand_char <- as.character(strand(tx_exons)[1])

  # Simple exon annotation (no CDS/UTR info for novel isoforms)
  exon_names <- paste0(start(tx_exons), "_", end(tx_exons), "_", strand_char, ":exon")

  # Junctions
  if (length(tx_exons) > 1) {
    junctions <- paste0(end(tx_exons[-length(tx_exons)]), "_", start(tx_exons[-1]), "_", strand_char)
    junctions_str <- paste(junctions, collapse = ";")
  } else {
    junctions_str <- NA
  }

  return(list(
    exons = paste(exon_names, collapse = ";"),
    junctions = junctions_str,
    n_exons = length(tx_exons)
  ))
}

get_cds_info <- function(gtf, transcript_id) {
  # Extract CDS length, protein length, and CDS genomic coordinates
  cds <- gtf[gtf$type == "CDS" & gtf$transcript_id == transcript_id]
  if (length(cds) == 0) {
    return(list(
      has_cds = FALSE,
      cds_length_nt = NA,
      cds_length_aa = NA,
      cds_start = NA,
      cds_end = NA
    ))
  }
  cds_length <- sum(width(cds))

  # CDS start and end are genomic coordinates (always start < end regardless of strand)
  cds_start <- min(start(cds))
  cds_end <- max(end(cds))

  return(list(
    has_cds = TRUE,
    cds_length_nt = cds_length,
    cds_length_aa = floor(cds_length / 3),
    cds_start = cds_start,
    cds_end = cds_end
  ))
}

calculate_utr_lengths <- function(tss, tes, cds_start, cds_end, strand) {
  # Calculate UTR lengths from TSS, TES, and CDS positions (strand-aware)
  # Returns NA if no CDS information available

  if (is.na(cds_start) || is.na(cds_end) || is.na(tss) || is.na(tes)) {
    return(list(utr5_length = NA, utr3_length = NA))
  }

  if (strand == "+") {
    # Plus strand: TSS < CDS_start < CDS_end < TES
    utr5_length <- max(0, cds_start - tss)
    utr3_length <- max(0, tes - cds_end)
  } else if (strand == "-") {
    # Minus strand: TES < CDS_start < CDS_end < TSS
    utr5_length <- max(0, tss - cds_end)
    utr3_length <- max(0, cds_start - tes)
  } else {
    return(list(utr5_length = NA, utr3_length = NA))
  }

  return(list(
    utr5_length = utr5_length,
    utr3_length = utr3_length
  ))
}

get_tss_tes <- function(gtf, transcript_id) {
  # Get transcription start and end sites
  tx <- gtf[gtf$type == "transcript" & gtf$transcript_id == transcript_id]
  if (length(tx) == 0) return(list(tss = NA, tes = NA))

  if (as.character(strand(tx)) == "+") {
    return(list(tss = start(tx), tes = end(tx)))
  } else {
    return(list(tss = end(tx), tes = start(tx)))
  }
}

# =============================================================================
# HELPER FUNCTIONS FOR PROTEIN TRANSLATION
# =============================================================================

get_cds_transcript_coords <- function(exons, cds_features, strand) {
  # Map genomic CDS boundaries to 1-based transcript-relative positions
  # exons: GRanges of exons for this transcript (ordered by genomic position)
  # cds_features: GRanges of CDS features for this transcript
  # strand: "+" or "-"
  # Returns list(tx_start, tx_end) — 1-based positions within the mRNA sequence

  if (length(cds_features) == 0) return(NULL)

  exons <- exons[order(start(exons))]

  # For minus strand, transcript order is reverse of genomic order
  if (strand == "-") {
    exon_order <- rev(seq_along(exons))
  } else {
    exon_order <- seq_along(exons)
  }

  # Build cumulative exon lengths in transcript order
  exon_widths <- width(exons)[exon_order]
  cumulative <- c(0, cumsum(exon_widths))

  # CDS genomic boundaries
  cds_start_genomic <- min(start(cds_features))
  cds_end_genomic <- max(end(cds_features))

  # Find transcript-relative CDS start
  tx_cds_start <- NA
  tx_cds_end <- NA

  for (idx in seq_along(exon_order)) {
    ei <- exon_order[idx]
    e_start <- start(exons)[ei]
    e_end <- end(exons)[ei]

    if (strand == "+") {
      # CDS start: first base of CDS
      if (is.na(tx_cds_start) && cds_start_genomic >= e_start && cds_start_genomic <= e_end) {
        offset_in_exon <- cds_start_genomic - e_start
        tx_cds_start <- cumulative[idx] + offset_in_exon + 1
      }
      # CDS end: last base of CDS
      if (is.na(tx_cds_end) && cds_end_genomic >= e_start && cds_end_genomic <= e_end) {
        offset_in_exon <- cds_end_genomic - e_start
        tx_cds_end <- cumulative[idx] + offset_in_exon + 1
      }
    } else {
      # Minus strand: CDS end (highest genomic coord) maps to transcript start
      if (is.na(tx_cds_start) && cds_end_genomic >= e_start && cds_end_genomic <= e_end) {
        offset_in_exon <- e_end - cds_end_genomic
        tx_cds_start <- cumulative[idx] + offset_in_exon + 1
      }
      # CDS start (lowest genomic coord) maps to transcript end
      if (is.na(tx_cds_end) && cds_start_genomic >= e_start && cds_start_genomic <= e_end) {
        offset_in_exon <- e_end - cds_start_genomic
        tx_cds_end <- cumulative[idx] + offset_in_exon + 1
      }
    }
  }

  if (is.na(tx_cds_start) || is.na(tx_cds_end)) return(NULL)

  # Sanity check: mapped CDS length should equal sum of CDS feature widths
  expected_cds_len <- sum(width(cds_features))
  actual_cds_len <- tx_cds_end - tx_cds_start + 1
  if (actual_cds_len != expected_cds_len) {
    warning(paste("CDS length mismatch - expected:", expected_cds_len,
                  "actual:", actual_cds_len))
    return(NULL)
  }

  return(list(tx_start = tx_cds_start, tx_end = tx_cds_end))
}

translate_cds <- function(transcript_seq, tx_cds_start, tx_cds_end) {
  # Extract CDS subsequence from transcript DNAString and translate to protein
  # transcript_seq: DNAString of the full transcript
  # tx_cds_start, tx_cds_end: 1-based transcript-relative CDS positions
  # Returns: AAString of protein sequence, or NA on failure

  tryCatch({
    cds_seq <- Biostrings::subseq(transcript_seq, start = tx_cds_start, end = tx_cds_end)
    protein <- Biostrings::translate(cds_seq)
    return(protein)
  }, error = function(e) {
    return(NA)
  })
}

fetch_uniprot_canonical <- function(gene_name) {
  # Fetch canonical Uniprot protein sequence for a human gene via REST API
  # Returns list(header, sequence) or NULL on failure

  encoded_gene <- utils::URLencode(gene_name, reserved = TRUE)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:",
    encoded_gene,
    "+AND+organism_id:9606+AND+reviewed:true&format=fasta"
  )

  tryCatch({
    old_timeout <- getOption("timeout")
    options(timeout = 30)
    on.exit(options(timeout = old_timeout), add = TRUE)
    response <- readLines(url, warn = FALSE)

    if (length(response) == 0) {
      message(paste("    WARNING: No Uniprot entry found for", gene_name))
      return(NULL)
    }

    # Parse FASTA: find header lines
    header_idx <- grep("^>", response)
    if (length(header_idx) == 0) return(NULL)

    # Take first (canonical) entry
    header <- response[header_idx[1]]
    if (length(header_idx) > 1) {
      seq_lines <- response[(header_idx[1] + 1):(header_idx[2] - 1)]
    } else {
      seq_lines <- response[(header_idx[1] + 1):length(response)]
    }
    sequence <- paste(seq_lines, collapse = "")

    return(list(header = header, sequence = sequence))
  }, error = function(e) {
    message(paste("    WARNING: Failed to fetch Uniprot sequence:", e$message))
    return(NULL)
  })
}

read_uniprot_from_file <- function(file_path, gene_name) {
  # Read a local Uniprot FASTA and find the entry matching gene_name
  # Returns list(header, sequence) or NULL if not found

  tryCatch({
    lines <- readLines(file_path, warn = FALSE)
    header_idx <- grep("^>", lines)

    if (length(header_idx) == 0) return(NULL)

    # Search headers for gene name (Uniprot format: >sp|ACCESSION|GENE_HUMAN ...)
    # Also check GN=GENE_NAME in the header
    for (i in seq_along(header_idx)) {
      hdr <- lines[header_idx[i]]
      # Match gene name in sp|...|GENE_HUMAN or GN=GENE patterns
      gene_pattern <- paste0("\\b", gene_name, "[_ ]|GN=", gene_name, "\\b")
      if (grepl(gene_pattern, hdr, ignore.case = TRUE, perl = TRUE)) {
        # Extract sequence lines
        if (i < length(header_idx)) {
          seq_lines <- lines[(header_idx[i] + 1):(header_idx[i + 1] - 1)]
        } else {
          seq_lines <- lines[(header_idx[i] + 1):length(lines)]
        }
        sequence <- paste(seq_lines, collapse = "")
        return(list(header = hdr, sequence = sequence))
      }
    }

    message(paste("    WARNING: Gene", gene_name, "not found in Uniprot file"))
    return(NULL)
  }, error = function(e) {
    message(paste("    WARNING: Failed to read Uniprot file:", e$message))
    return(NULL)
  })
}

# =============================================================================
# VALIDATION CHECKS
# =============================================================================

message("====================================================================")
message("Gene Isoform Annotation Pipeline")
message("====================================================================")
message(paste("Gene:", GENE_INPUT))
message(paste("Date:", Sys.Date()))
message("")

# Check required files
message("Checking required files...")
files_to_check <- c(
  "GENCODE gene lookup" = GENCODE_GENE_LOOKUP,
  "GENCODE GTF (indexed)" = GENCODE_GTF_INDEXED,
  "GENCODE GTF index" = paste0(GENCODE_GTF_INDEXED, ".tbi"),
  "SQANTI GTF (indexed)" = SQANTI_GTF_INDEXED,
  "SQANTI GTF index" = paste0(SQANTI_GTF_INDEXED, ".tbi"),
  "SQANTI Classification" = SQANTI_CLASSIFICATION
)

if (INCLUDE_EXPRESSION) {
  files_to_check["DGEList RDS"] <- DGELIST_RDS
}

if (OUTPUT_FASTA || OUTPUT_PROTEIN) {
  files_to_check["GENCODE FASTA"] <- GENCODE_FASTA
  files_to_check["SQANTI FASTA"] <- SQANTI_FASTA
}

if (OUTPUT_PROTEIN && file.exists(SQANTI_PROTEIN_FASTA)) {
  files_to_check["SQANTI Protein FASTA"] <- SQANTI_PROTEIN_FASTA
}

missing_files <- c()
for (file_desc in names(files_to_check)) {
  file_path <- files_to_check[file_desc]
  if (file.exists(file_path)) {
    message(paste("  [OK]", file_desc))
  } else {
    message(paste("  [MISSING]", file_desc, "-", file_path))
    missing_files <- c(missing_files, file_path)
  }
}

if (length(missing_files) > 0) {
  stop(paste("ERROR: Missing required files. Please run preprocess_sqanti_tabix.sh first to create indexed GTF files."))
}

# Check for tabix command
if (Sys.which("tabix") == "") {
  stop("ERROR: tabix command not found. Install with: conda install -c bioconda htslib")
}

# Check for Biostrings if FASTA or protein output requested
if (OUTPUT_FASTA || OUTPUT_PROTEIN) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("ERROR: Biostrings package required for FASTA/protein output.\nInstall with: BiocManager::install('Biostrings')")
  }
}

# Validate --uniprot_file if provided
if (!is.null(UNIPROT_FILE) && !file.exists(UNIPROT_FILE)) {
  stop(paste("ERROR: Uniprot file not found:", UNIPROT_FILE))
}

message("")
message("All checks passed!")
message("")

# =============================================================================
# STEP 1: FIND GENE ID
# =============================================================================

message("====================================================================")
message("STEP 1: Finding gene in GENCODE")
message("====================================================================")

if (SEARCH_BY == "name") {
  message(paste("  Searching by gene name:", GENE_NAME))
  gene_lines <- system(
    paste0("grep 'gene_name \"", GENE_NAME, "\"' ", GENCODE_GENE_LOOKUP,
           " | awk -F'\\t' '$3 == \"gene\"'"),
    intern = TRUE
  )

  if (length(gene_lines) == 0) {
    stop(paste("ERROR: Gene name", GENE_NAME, "not found in GENCODE GTF"))
  }

  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines[1])
  gene_name_found <- sub('.*gene_name "([^"]+)".*', '\\1', gene_lines[1])
  GENE_NAME <- gene_name_found  # Use canonical GENCODE name for output files

  message(paste("  Found gene ID:", gene_id))
  message(paste("  Gene name:", gene_name_found))

} else {
  message(paste("  Searching by Ensembl ID:", GENE_ID))
  gene_lines <- system(
    paste0("grep 'gene_id \"", GENE_ID, "' ", GENCODE_GENE_LOOKUP, " | grep -w gene | head -1"),
    intern = TRUE
  )

  if (length(gene_lines) == 0) {
    stop(paste("ERROR: Gene ID", GENE_ID, "not found in GENCODE GTF"))
  }

  gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_lines[1])
  gene_name_found <- sub('.*gene_name "([^"]+)".*', '\\1', gene_lines[1])

  GENE_NAME <- gene_name_found  # Set for output file naming

  message(paste("  Found gene ID:", gene_id))
  message(paste("  Gene name:", gene_name_found))
}

# Parse gene coordinates from the gene line for tabix query
gene_parts <- strsplit(gene_lines[1], "\t")[[1]]
gene_chr <- gene_parts[1]
gene_start <- as.numeric(gene_parts[4])
gene_end <- as.numeric(gene_parts[5])

message(paste("  Gene coordinates:", paste0(gene_chr, ":", gene_start, "-", gene_end)))

# Set output filenames now that GENE_NAME is determined
OUTPUT_TSV <- paste0("isoform_annotation_", GENE_NAME, ".tsv")
OUTPUT_GTF_FILE <- paste0("isoform_annotation_", GENE_NAME, ".gtf")
OUTPUT_FASTA_FILE <- paste0("isoform_annotation_", GENE_NAME, ".fasta")
OUTPUT_PROTEIN_FILE <- paste0("protein_sequences_", GENE_NAME, ".fasta")

message("")

# =============================================================================
# STEP 2: EXTRACT GENCODE ISOFORMS
# =============================================================================

message("====================================================================")
message("STEP 2: Extracting GENCODE isoforms")
message("====================================================================")

# Query tabix-indexed GENCODE GTF by coordinates (fast: <1 second)
message("  Querying tabix-indexed GENCODE GTF...")
tabix_cmd <- sprintf("tabix %s %s:%d-%d", GENCODE_GTF_INDEXED, gene_chr, gene_start, gene_end)
gencode_lines <- system(tabix_cmd, intern = TRUE)

if (length(gencode_lines) == 0) {
  stop(paste("ERROR: No features found for gene coordinates:", paste0(gene_chr, ":", gene_start, "-", gene_end)))
}

# Filter for this specific gene (tabix returns all features in region, may include overlapping genes)
gencode_lines <- gencode_lines[grepl(paste0('gene_id "', gene_id), gencode_lines)]

message(paste("  Retrieved", length(gencode_lines), "features from indexed GTF"))

# Parse into GRanges
temp_gtf <- tempfile(fileext = ".gtf")
writeLines(gencode_lines, temp_gtf)
gencode <- import(temp_gtf)
unlink(temp_gtf)

gencode_transcripts <- gencode[gencode$type == "transcript" & gencode$gene_id == gene_id]
message(paste("  Found", length(gencode_transcripts), "GENCODE transcripts"))

gencode_result_list <- vector("list", length(gencode_transcripts))

message("  Extracting detailed annotations...")
for (i in seq_along(gencode_transcripts)) {
  tx <- gencode_transcripts[i]
  tx_id <- tx$transcript_id

  exon_info <- get_exon_structure_annotated(gencode, tx_id)
  cds_info <- get_cds_info(gencode, tx_id)
  tss_tes <- get_tss_tes(gencode, tx_id)

  # Calculate UTR lengths from TSS/TES/CDS positions (strand-aware)
  utr_info <- calculate_utr_lengths(
    tss_tes$tss,
    tss_tes$tes,
    cds_info$cds_start,
    cds_info$cds_end,
    as.character(strand(tx))
  )

  is_protein_coding <- FALSE
  coding_source <- NA
  if (!is.na(tx$transcript_type) && tx$transcript_type == "protein_coding") {
    is_protein_coding <- TRUE
    coding_source <- "GENCODE_biotype"
  } else if (cds_info$has_cds) {
    is_protein_coding <- TRUE
    coding_source <- "GENCODE_CDS"
  }

  gencode_result_list[[i]] <- data.frame(
    isoform_id = tx_id,
    isoform_name = ifelse(!is.na(tx$transcript_name), tx$transcript_name, tx_id),
    source = "GENCODE",
    chromosome = as.character(seqnames(tx)),
    strand = as.character(strand(tx)),
    tss = tss_tes$tss,
    tes = tss_tes$tes,
    transcript_length = width(tx),
    n_exons = exon_info$n_exons,
    exons = exon_info$exons,
    junctions = exon_info$junctions,
    cds_start = cds_info$cds_start,
    cds_end = cds_info$cds_end,
    cds_length_nt = cds_info$cds_length_nt,
    cds_length_aa = cds_info$cds_length_aa,
    in_frame = if (cds_info$has_cds) as.integer(cds_info$cds_length_nt %% 3 == 0) else NA_integer_,
    utr5_length = utr_info$utr5_length,
    utr3_length = utr_info$utr3_length,
    is_protein_coding = is_protein_coding,
    coding_source = coding_source,
    stringsAsFactors = FALSE
  )
}
gencode_results <- bind_rows(gencode_result_list)

message(paste("  Extracted", nrow(gencode_results), "GENCODE isoforms with full annotations"))
message("")

# =============================================================================
# STEP 3: EXTRACT SQANTI ISOFORMS
# =============================================================================

message("====================================================================")
message("STEP 3: Extracting SQANTI isoforms")
message("====================================================================")

# Part A: Get isoform metadata from SQANTI classification
message("  Loading SQANTI classification...")
temp_sqanti <- tempfile(fileext = ".txt")
system(paste0("head -1 ", SQANTI_CLASSIFICATION, " > ", temp_sqanti))
system(paste0("grep '", gene_id, "' ", SQANTI_CLASSIFICATION, " >> ", temp_sqanti))

sqanti_class <- read_tsv(temp_sqanti, show_col_types = FALSE)
unlink(temp_sqanti)

sqanti_target <- sqanti_class %>%
  filter(grepl(paste0("\\b", gene_id, "\\b"), associated_gene, ignore.case = TRUE))

message(paste("  Found", nrow(sqanti_target), "SQANTI isoforms"))

sqanti_results <- data.frame()
sqanti_features <- NULL

if (nrow(sqanti_target) > 0) {
  # Part B: Get exon structures from tabix-indexed SQANTI GTF
  message("  Querying tabix-indexed SQANTI GTF for exon structures...")

  # Query tabix by gene coordinates
  tabix_cmd <- sprintf("tabix %s %s:%d-%d", SQANTI_GTF_INDEXED, gene_chr, gene_start, gene_end)
  sqanti_lines <- system(tabix_cmd, intern = TRUE)

  if (length(sqanti_lines) > 0) {
    message(paste("    Retrieved", length(sqanti_lines), "GTF lines"))

    # Parse tabix output into GRanges
    temp_gtf <- tempfile(fileext = ".gtf")
    writeLines(sqanti_lines, temp_gtf)
    sqanti_features <- import(temp_gtf)
    unlink(temp_gtf)

    message("    Parsed SQANTI GTF features")
  } else {
    message("    No SQANTI features found in gene region")
  }

  # Process each SQANTI isoform
  message("  Extracting SQANTI isoform annotations...")
  sqanti_result_list <- vector("list", nrow(sqanti_target))
  for (i in seq_len(nrow(sqanti_target))) {
    row <- sqanti_target[i, ]
    pb_id <- row$isoform

    # Get exon structure and CDS info from GTF if available
    exon_info <- list(exons = NA, junctions = NA, n_exons = row$exons)
    tss_tes <- list(tss = NA, tes = NA)
    cds_info <- list(has_cds = FALSE, cds_start = NA, cds_end = NA, cds_length_nt = NA, cds_length_aa = NA)

    if (!is.null(sqanti_features)) {
      # Check if this isoform is in the GTF
      tx_features <- sqanti_features[sqanti_features$transcript_id == pb_id]

      if (length(tx_features) > 0) {
        exon_info <- get_exon_structure_simple(sqanti_features, pb_id)
        tss_tes <- get_tss_tes(sqanti_features, pb_id)
        cds_info <- get_cds_info(sqanti_features, pb_id)
      }
    }

    # Use SQANTI classification for CDS/protein if GTF doesn't have it
    if (!cds_info$has_cds) {
      cds_info$cds_length_nt <- ifelse(!is.na(row$CDS_length), row$CDS_length, NA)
      cds_info$cds_length_aa <- ifelse(!is.na(row$ORF_length), row$ORF_length, NA)
    }

    # Calculate UTR lengths from TSS/TES/CDS if available
    utr_info <- calculate_utr_lengths(
      tss_tes$tss,
      tss_tes$tes,
      cds_info$cds_start,
      cds_info$cds_end,
      row$strand
    )

    # Determine coding status
    is_protein_coding <- FALSE
    coding_source <- NA
    if (!is.na(row$coding) && row$coding == "coding") {
      is_protein_coding <- TRUE
      coding_source <- "SQANTI"
    }

    sqanti_result_list[[i]] <- data.frame(
      isoform_id = pb_id,
      isoform_name = pb_id,
      source = "SQANTI",
      chromosome = row$chrom,
      strand = row$strand,
      tss = tss_tes$tss,
      tes = tss_tes$tes,
      transcript_length = row$length,
      n_exons = exon_info$n_exons,
      exons = exon_info$exons,
      junctions = exon_info$junctions,
      cds_start = cds_info$cds_start,
      cds_end = cds_info$cds_end,
      cds_length_nt = cds_info$cds_length_nt,
      cds_length_aa = cds_info$cds_length_aa,
      in_frame = if (!is.na(cds_info$cds_length_nt)) as.integer(cds_info$cds_length_nt %% 3 == 0) else NA_integer_,
      utr5_length = utr_info$utr5_length,
      utr3_length = utr_info$utr3_length,
      is_protein_coding = is_protein_coding,
      coding_source = coding_source,
      stringsAsFactors = FALSE
    )
  }
  sqanti_results <- bind_rows(sqanti_result_list)

  message(paste("  Extracted", nrow(sqanti_results), "SQANTI isoforms"))
}

message("")

# =============================================================================
# STEP 4: CALCULATE EXPRESSION (OPTIONAL)
# =============================================================================

# Combine GENCODE and SQANTI results
all_isoforms <- rbind(gencode_results, sqanti_results)

if (INCLUDE_EXPRESSION) {
  message("====================================================================")
  message("STEP 4: Calculating expression and read counts")
  message("====================================================================")
  message(paste("  Total isoforms:", nrow(all_isoforms)))

  # Load DGEList
  message("  Loading DGEList...")
  dge_isoform <- readRDS(DGELIST_RDS)

  # Get normalized CPM and raw counts
  message("  Calculating normalized CPM...")
  cpm_normalized <- cpm(dge_isoform, normalized = TRUE, log = FALSE)
  raw_counts <- dge_isoform$counts

  # Get transcript IDs from genes component
  if (!TRANSCRIPT_ID_COLUMN %in% colnames(dge_isoform$genes)) {
    stop("ERROR: '", TRANSCRIPT_ID_COLUMN, "' column not found in DGEList$genes. ",
         "Available columns: ", paste(colnames(dge_isoform$genes), collapse = ", "),
         "\nSet TRANSCRIPT_ID_COLUMN in config.R to the correct column name.")
  }
  transcript_ids <- dge_isoform$genes[[TRANSCRIPT_ID_COLUMN]]
  rownames(cpm_normalized) <- transcript_ids
  rownames(raw_counts) <- transcript_ids

  # Use rownames as sample IDs (standard DGEList convention)
  sample_info <- dge_isoform$samples
  sample_ids <- rownames(sample_info)

  # Build sample groups from STRATIFY_BY configuration
  if (is.null(STRATIFY_BY)) {
    # No stratification: single group across all samples
    group_labels <- "mean"
    sample_group <- setNames(rep("mean", length(sample_ids)), sample_ids)
    message("  Calculating expression and read counts (unstratified)...")
  } else {
    # Validate STRATIFY_BY columns exist in sample metadata
    missing_cols <- setdiff(STRATIFY_BY, colnames(sample_info))
    if (length(missing_cols) > 0) {
      stop("ERROR: STRATIFY_BY column(s) not found in DGEList$samples: ",
           paste(missing_cols, collapse = ", "),
           "\nAvailable columns: ", paste(colnames(sample_info), collapse = ", "))
    }

    # Get ordered levels per column, applying LEVEL_LABELS where provided
    level_orders <- list()
    for (col in STRATIFY_BY) {
      raw_levels <- if (is.factor(sample_info[[col]])) {
        levels(sample_info[[col]])
      } else {
        unique(as.character(sample_info[[col]]))
      }

      if (!is.null(LEVEL_LABELS) && col %in% names(LEVEL_LABELS) &&
          !is.null(LEVEL_LABELS[[col]])) {
        label_map <- LEVEL_LABELS[[col]]
        unmapped <- setdiff(raw_levels, names(label_map))
        if (length(unmapped) > 0) {
          stop("ERROR: LEVEL_LABELS for '", col, "' is missing mappings for: ",
               paste(unmapped, collapse = ", "),
               "\nAll levels must be mapped when LEVEL_LABELS is provided for a column.")
        }
        level_orders[[col]] <- unname(label_map[names(label_map) %in% raw_levels])
      } else {
        level_orders[[col]] <- raw_levels
      }
    }

    # Build group label for each sample (vectorized)
    label_parts <- lapply(STRATIFY_BY, function(col) {
      vals <- as.character(sample_info[[col]])
      if (!is.null(LEVEL_LABELS) && col %in% names(LEVEL_LABELS) &&
          !is.null(LEVEL_LABELS[[col]])) {
        vals <- LEVEL_LABELS[[col]][vals]
      }
      vals
    })
    sample_group <- setNames(
      do.call(paste, c(label_parts, list(sep = "_"))),
      sample_ids
    )

    # Iterate groups in deterministic order (first STRATIFY_BY column varies fastest)
    group_grid <- do.call(expand.grid, c(level_orders, list(stringsAsFactors = FALSE)))
    group_labels <- apply(group_grid, 1, paste, collapse = "_")

    message("  Calculating expression and read counts per group...")
  }

  # Calculate mean CPM and total raw read counts per group
  for (group_label in group_labels) {
    group_samples <- sample_ids[sample_group == group_label]

    expr_col  <- paste0("expr_", group_label)
    reads_col <- paste0("total_reads_", group_label)

    if (length(group_samples) > 0) {
      # Mean TMM-normalized CPM across replicates
      cpm_means <- rowMeans(cpm_normalized[, group_samples, drop = FALSE])
      all_isoforms[[expr_col]] <- cpm_means[all_isoforms$isoform_id]

      # Sum of raw counts across replicates
      read_sums <- rowSums(raw_counts[, group_samples, drop = FALSE])
      all_isoforms[[reads_col]] <- read_sums[all_isoforms$isoform_id]

      n_expressed <- sum(!is.na(all_isoforms[[expr_col]]) & all_isoforms[[expr_col]] > 0)
      n_missing   <- sum(is.na(all_isoforms[[expr_col]]))
      message(paste("    ", group_label, ":", n_expressed, "isoforms with expression > 0"))
      if (n_missing > 0) {
        message(paste("      Note:", n_missing, "isoforms not found in DGEList"))
      }
    } else {
      all_isoforms[[expr_col]]  <- NA
      all_isoforms[[reads_col]] <- NA
      message(paste("    ", group_label, ": No samples found"))
    }
  }

  message("")
} else {
  message("====================================================================")
  message("STEP 4: Skipping expression calculation (INCLUDE_EXPRESSION = FALSE)")
  message("====================================================================")
  message(paste("  Total isoforms:", nrow(all_isoforms)))
  message("")
}

# =============================================================================
# STEP 5: WRITE OUTPUT FILES
# =============================================================================

message("====================================================================")
message("STEP 5: Writing output files")
message("====================================================================")

# Round CPM expression values to 2 significant digits (if included)
if (INCLUDE_EXPRESSION) {
  expr_cols <- grep("^expr_", colnames(all_isoforms), value = TRUE)
  for (col in expr_cols) {
    all_isoforms[[col]] <- signif(all_isoforms[[col]], 2)
  }
}

# Replace NA with 0 in expression and read count columns
if (INCLUDE_EXPRESSION) {
  all_expr_reads_cols <- grep("^(expr_|total_reads_)", colnames(all_isoforms), value = TRUE)
  for (col in all_expr_reads_cols) {
    all_isoforms[[col]][is.na(all_isoforms[[col]])] <- 0
  }
}

# Add total reads across all samples and conditions
if (INCLUDE_EXPRESSION) {
  reads_cols <- grep("^total_reads_", colnames(all_isoforms), value = TRUE)
  all_isoforms$total_reads_all_samples <- rowSums(
    all_isoforms[, reads_cols, drop = FALSE]
  )
}

# Optional: Write GTF
if (OUTPUT_GTF) {
  message(paste("  Writing GTF:", OUTPUT_GTF_FILE))

  # Combine GENCODE and SQANTI GTF features
  output_gtf <- gencode

  if (!is.null(sqanti_features) && length(sqanti_features) > 0) {
    # Filter to only isoforms belonging to this gene (exclude overlapping genes)
    sqanti_gene_ids <- sqanti_target$isoform
    sqanti_features_filtered <- sqanti_features[
      sqanti_features$transcript_id %in% sqanti_gene_ids
    ]
    output_gtf <- c(output_gtf, sqanti_features_filtered)
  }

  # Add chr prefix to seqnames for UCSC browser compatibility
  seqlvls <- seqlevels(output_gtf)
  seqlevels(output_gtf) <- ifelse(startsWith(seqlvls, "chr"), seqlvls, paste0("chr", seqlvls))

  export(output_gtf, OUTPUT_GTF_FILE, format = "gtf")
  message(paste("    Wrote", length(output_gtf), "GTF features"))
}

# Optional: Write FASTA
if (OUTPUT_FASTA) {
  message(paste("  Writing FASTA:", OUTPUT_FASTA_FILE))
  message("    (This may take 1-2 minutes for large FASTA files)")

  library(Biostrings)

  # Separate GENCODE and SQANTI isoforms
  gencode_ids <- all_isoforms$isoform_id[all_isoforms$source == "GENCODE"]
  sqanti_ids <- all_isoforms$isoform_id[all_isoforms$source == "SQANTI"]

  # Read GENCODE FASTA and extract GENCODE isoforms
  message("    Reading GENCODE transcripts FASTA...")
  gencode_fasta <- readDNAStringSet(GENCODE_FASTA)

  # GENCODE FASTA has pipe-delimited headers like: ENST00000832824.1|ENSG...|...
  # Extract transcript ID (first field before pipe)
  gencode_fasta_ids <- sub("\\|.*", "", names(gencode_fasta))

  # Match by transcript ID
  matches_idx <- gencode_fasta_ids %in% gencode_ids
  gencode_matches <- gencode_fasta[matches_idx]

  # Simplify names to just transcript ID
  names(gencode_matches) <- gencode_fasta_ids[matches_idx]

  message(paste("      Found", length(gencode_matches), "GENCODE sequences"))

  # Read SQANTI FASTA and extract SQANTI isoforms
  message("    Reading SQANTI transcripts FASTA...")
  sqanti_fasta <- readDNAStringSet(SQANTI_FASTA)
  sqanti_matches <- sqanti_fasta[names(sqanti_fasta) %in% sqanti_ids]
  message(paste("      Found", length(sqanti_matches), "SQANTI sequences"))

  # Combine sequences
  all_sequences <- c(gencode_matches, sqanti_matches)

  if (length(all_sequences) > 0) {
    writeXStringSet(all_sequences, OUTPUT_FASTA_FILE, format = "fasta")
    message(paste("    Wrote", length(all_sequences), "total sequences to", OUTPUT_FASTA_FILE))
  } else {
    message("    WARNING: No matching sequences found in FASTA files")
  }
}

# Optional: Write protein FASTA
if (OUTPUT_PROTEIN) {
  message(paste("  Writing protein FASTA:", OUTPUT_PROTEIN_FILE))

  library(Biostrings)

  protein_sequences <- AAStringSet()

  # Step 1: Get canonical Uniprot sequence
  message("    Fetching canonical Uniprot sequence...")
  if (!is.null(UNIPROT_FILE)) {
    uniprot <- read_uniprot_from_file(UNIPROT_FILE, GENE_NAME)
  } else {
    uniprot <- fetch_uniprot_canonical(GENE_NAME)
  }

  if (!is.null(uniprot)) {
    canonical_aa <- AAStringSet(uniprot$sequence)
    names(canonical_aa) <- sub("^>", "", uniprot$header)
    protein_sequences <- c(protein_sequences, canonical_aa)
    message(paste("      Uniprot canonical:", nchar(uniprot$sequence), "aa"))
  } else {
    message("      WARNING: No Uniprot canonical sequence available")
  }

  # Step 2: Load transcript sequences (reuse from FASTA step if available, otherwise load fresh)
  if (!exists("gencode_fasta") || !exists("sqanti_fasta")) {
    message("    Loading transcript FASTA files for translation...")
    gencode_fasta <- readDNAStringSet(GENCODE_FASTA)
    gencode_fasta_ids <- sub("\\|.*", "", names(gencode_fasta))
    sqanti_fasta <- readDNAStringSet(SQANTI_FASTA)
  }

  # Pre-load SQANTI protein FASTA for fallback translation
  sqanti_protein_fasta <- NULL
  if (file.exists(SQANTI_PROTEIN_FASTA)) {
    tryCatch({
      message("    Loading SQANTI protein FASTA for fallback...")
      sqanti_protein_fasta <- readAAStringSet(SQANTI_PROTEIN_FASTA)
      # Headers may be tab-delimited (e.g., "PB.14156.12\tPB.14156.12.p2|...") — keep only first field
      names(sqanti_protein_fasta) <- sub("\t.*", "", names(sqanti_protein_fasta))
      message(paste("      Loaded", length(sqanti_protein_fasta), "protein sequences"))
    }, error = function(e) {
      message(paste("    WARNING: Failed to read SQANTI protein FASTA:", e$message))
    })
  }

  # Step 3: Translate protein-coding isoforms
  protein_coding <- all_isoforms[all_isoforms$is_protein_coding == TRUE, ]
  message(paste("    Translating", nrow(protein_coding), "protein-coding isoforms..."))

  n_translated <- 0
  n_failed <- 0
  translated_proteins <- list()  # Store translations for UniProt comparison

  for (i in seq_len(nrow(protein_coding))) {
    iso <- protein_coding[i, ]
    iso_id <- iso$isoform_id
    iso_source <- iso$source

    translated <- NULL

    if (iso_source == "GENCODE") {
      # Get transcript sequence from GENCODE FASTA
      fasta_match <- which(gencode_fasta_ids == iso_id)
      if (length(fasta_match) == 0) {
        n_failed <- n_failed + 1
        next
      }
      tx_seq <- gencode_fasta[[fasta_match[1]]]

      # Get CDS features from GTF
      cds_feats <- gencode[gencode$type == "CDS" & gencode$transcript_id == iso_id]
      tx_exons <- gencode[gencode$type == "exon" & gencode$transcript_id == iso_id]
      if (length(cds_feats) == 0 || length(tx_exons) == 0) {
        n_failed <- n_failed + 1
        next
      }

      strand_char <- as.character(strand(tx_exons)[1])
      coords <- get_cds_transcript_coords(tx_exons, cds_feats, strand_char)
      if (is.null(coords)) {
        n_failed <- n_failed + 1
        next
      }

      translated <- translate_cds(tx_seq, coords$tx_start, coords$tx_end)

    } else if (iso_source == "SQANTI") {
      # Try SQANTI GTF CDS features first
      has_sqanti_cds <- FALSE
      if (!is.null(sqanti_features)) {
        cds_feats <- sqanti_features[sqanti_features$type == "CDS" & sqanti_features$transcript_id == iso_id]
        tx_exons <- sqanti_features[sqanti_features$type == "exon" & sqanti_features$transcript_id == iso_id]
        if (length(cds_feats) > 0 && length(tx_exons) > 0) {
          has_sqanti_cds <- TRUE
          # Get transcript sequence
          sq_match <- which(names(sqanti_fasta) == iso_id)
          if (length(sq_match) > 0) {
            tx_seq <- sqanti_fasta[[sq_match[1]]]
            strand_char <- as.character(strand(tx_exons)[1])
            coords <- get_cds_transcript_coords(tx_exons, cds_feats, strand_char)
            if (!is.null(coords)) {
              translated <- translate_cds(tx_seq, coords$tx_start, coords$tx_end)
            }
          }
        }
      }

      # Fallback: try pre-loaded sqanti3_corrected.faa
      if (is.null(translated) && !is.null(sqanti_protein_fasta)) {
        faa_match <- which(names(sqanti_protein_fasta) == iso_id)
        if (length(faa_match) > 0) {
          translated <- sqanti_protein_fasta[[faa_match[1]]]
        }
      }
    }

    if (!is.null(translated) && !identical(translated, NA)) {
      aa_entry <- AAStringSet(translated)
      names(aa_entry) <- iso_id
      protein_sequences <- c(protein_sequences, aa_entry)
      translated_proteins[[iso_id]] <- as.character(translated)
      n_translated <- n_translated + 1
    } else {
      n_failed <- n_failed + 1
    }
  }

  message(paste("      Translated:", n_translated, "isoforms"))
  if (n_failed > 0) {
    message(paste("      Failed:", n_failed, "isoforms"))
  }

  # Step 4: Compare translated proteins to UniProt canonical
  all_isoforms$frame_matches_uniprot <- NA
  all_isoforms$match_pct_uniprot <- NA_real_
  all_isoforms$matched_range_uniprot <- NA_character_
  if (!is.null(uniprot) && length(translated_proteins) > 0) {
    uniprot_seq <- uniprot$sequence
    uniprot_len <- nchar(uniprot_seq)
    message("    Comparing translated proteins to UniProt canonical...")
    n_match <- 0
    n_no_match <- 0
    win_size <- 20

    for (iso_id in names(translated_proteins)) {
      # Strip trailing stop codon (*) for comparison — UniProt doesn't include it
      iso_seq <- sub("\\*$", "", translated_proteins[[iso_id]])
      is_match <- (iso_seq == uniprot_seq) || grepl(iso_seq, uniprot_seq, fixed = TRUE)
      row_idx <- which(all_isoforms$isoform_id == iso_id)

      if (length(row_idx) == 0) next

      all_isoforms$frame_matches_uniprot[row_idx] <- is_match
      if (is_match) n_match <- n_match + 1 else n_no_match <- n_no_match + 1

      # Window-based sequence comparison
      iso_chars <- strsplit(iso_seq, "")[[1]]
      iso_len <- length(iso_chars)

      if (iso_len < win_size) {
        # Too short for window analysis — use simple substring check
        all_isoforms$match_pct_uniprot[row_idx] <- if (grepl(iso_seq, uniprot_seq, fixed = TRUE)) 100 else 0
        if (grepl(iso_seq, uniprot_seq, fixed = TRUE)) {
          uni_pos <- regexpr(iso_seq, uniprot_seq, fixed = TRUE)
          all_isoforms$matched_range_uniprot[row_idx] <- paste0(uni_pos, "-", uni_pos + iso_len - 1)
        }
        next
      }

      n_windows <- iso_len - win_size + 1
      matches <- logical(n_windows)
      for (w in seq_len(n_windows)) {
        window <- paste0(iso_chars[w:(w + win_size - 1)], collapse = "")
        matches[w] <- grepl(window, uniprot_seq, fixed = TRUE)
      }

      pct <- round(100 * sum(matches) / n_windows, 1)
      all_isoforms$match_pct_uniprot[row_idx] <- pct

      # Find longest contiguous matching stretch and its UniProt position
      if (any(matches)) {
        rle_m <- rle(matches)
        longest_idx <- which.max(rle_m$lengths * rle_m$values)
        run_start <- sum(rle_m$lengths[seq_len(longest_idx - 1)]) + 1
        run_len <- rle_m$lengths[longest_idx]
        # Contiguous stretch spans from run_start to run_start + run_len - 1 + win_size - 1
        stretch_end <- run_start + run_len - 1 + win_size - 1
        stretch_seq <- paste0(iso_chars[run_start:stretch_end], collapse = "")
        uni_pos <- regexpr(stretch_seq, uniprot_seq, fixed = TRUE)
        if (uni_pos > 0) {
          all_isoforms$matched_range_uniprot[row_idx] <- paste0(uni_pos, "-", uni_pos + nchar(stretch_seq) - 1)
        }
      }
    }
    message(paste("      Exact/subset match:", n_match, "| Partial:", n_no_match,
                  "| Not translated:", sum(is.na(all_isoforms$frame_matches_uniprot) & all_isoforms$is_protein_coding)))
  }

  if (length(protein_sequences) > 0) {
    writeXStringSet(protein_sequences, OUTPUT_PROTEIN_FILE, format = "fasta")
    message(paste("    Wrote", length(protein_sequences), "protein sequences to", OUTPUT_PROTEIN_FILE))
  } else {
    message("    WARNING: No protein sequences generated")
  }
}

# Write TSV (after all columns including frame_matches_uniprot are populated)
message(paste("  Writing TSV:", OUTPUT_TSV))
write_tsv(all_isoforms, OUTPUT_TSV)
message(paste("    Wrote", nrow(all_isoforms), "isoforms"))

message("")

# =============================================================================
# SUMMARY
# =============================================================================

message("====================================================================")
message("SUMMARY")
message("====================================================================")
message(paste("Gene:", GENE_NAME, "(", gene_id, ")"))
message("")
message("Isoform counts:")
message(paste("  GENCODE:", sum(all_isoforms$source == "GENCODE")))
message(paste("  SQANTI:", sum(all_isoforms$source == "SQANTI")))
message(paste("  Total:", nrow(all_isoforms)))
message("")
message(paste("  Protein coding:", sum(all_isoforms$is_protein_coding, na.rm = TRUE)))
message(paste("  Non-coding:", sum(!all_isoforms$is_protein_coding, na.rm = TRUE)))
message("")
message("Output files:")
message(paste("  ", OUTPUT_TSV))
if (OUTPUT_GTF) message(paste("  ", OUTPUT_GTF_FILE))
if (OUTPUT_FASTA) message(paste("  ", OUTPUT_FASTA_FILE))
if (OUTPUT_PROTEIN) message(paste("  ", OUTPUT_PROTEIN_FILE))
message("")
message("====================================================================")
message("DONE!")
message("====================================================================")
