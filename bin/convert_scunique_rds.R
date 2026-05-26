#!/usr/bin/env Rscript
# convert_scunique_rds.R
#
# Converts an scUnique results directory (RDS files) into the plain-text
# formats expected by the sc_clone_mutations pipeline:
#   - MEDICC2 tree → Newick (.new)
#   - Unique events → TSV (cell_id, chr, start, end, event_type, copy_number)
#   - Cell metadata → CSV (cell_id, sample_id, patient_id, ploidy, pass_qc)
#
# Usage:
#   Rscript convert_scunique_rds.R <scunique_results_dir> <output_dir> [<sample_id>]

if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape", repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript convert_scunique_rds.R <scunique_dir> <output_dir> [<sample_id>] [<bam_dir>]")
}

input_dir  <- args[1]
output_dir <- args[2]

# Auto-detect sample prefix from RDS files
rds_files <- list.files(input_dir, pattern = "\\.RDS$", full.names = FALSE)
prefix <- unique(sub("\\.[^.]+\\.RDS$", "", rds_files))[1]
sample_id <- if (length(args) >= 3 && args[3] != "") args[3] else prefix
bam_dir   <- if (length(args) >= 4 && args[4] != "") args[4] else NULL

cat("Input dir :", input_dir, "\n")
cat("Prefix    :", prefix, "\n")
cat("Sample ID :", sample_id, "\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ─── Helper: find RDS file by pattern ─────────────────────────────────────────
find_rds <- function(pattern) {
  f <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) return(NULL)
  f[1]
}

# ─── 1. MEDICC2 tree → Newick ─────────────────────────────────────────────────
tree_file <- find_rds("medicc_tree\\.RDS$")
if (is.null(tree_file)) tree_file <- find_rds("tree\\.RDS$")

if (!is.null(tree_file)) {
  cat("Loading tree:", tree_file, "\n")
  tree_obj <- readRDS(tree_file)
  
  if (inherits(tree_obj, "phylo")) {
    out_tree <- file.path(output_dir, "medicc2_tree.new")
    write.tree(tree_obj, file = out_tree)
    cat("  -> Written:", out_tree, "\n")
  } else if (inherits(tree_obj, "hclust") || inherits(tree_obj, "dendrogram")) {
    tree_phylo <- as.phylo(as.hclust(tree_obj))
    out_tree <- file.path(output_dir, "medicc2_tree.new")
    write.tree(tree_phylo, file = out_tree)
    cat("  -> Converted hclust/dendrogram to phylo. Written:", out_tree, "\n")
  } else {
    cat("  WARNING: tree is class", class(tree_obj), "- attempting as.phylo\n")
    tryCatch({
      tree_phylo <- as.phylo(tree_obj)
      out_tree <- file.path(output_dir, "medicc2_tree.new")
      write.tree(tree_phylo, file = out_tree)
      cat("  -> Written:", out_tree, "\n")
    }, error = function(e) {
      cat("  ERROR: Cannot convert tree:", e$message, "\n")
    })
  }
} else {
  cat("WARNING: No medicc_tree.RDS or tree.RDS found\n")
}

# ─── 2. Unique events → TSV ──────────────────────────────────────────────────
events_file <- find_rds("uniqueEvents\\.RDS$")
if (is.null(events_file)) events_file <- find_rds("df_events\\.RDS$")

if (!is.null(events_file)) {
  cat("Loading events:", events_file, "\n")
  events_obj <- readRDS(events_file)
  
  if (is.data.frame(events_obj)) {
    events_df <- events_obj
  } else if (is.list(events_obj) && !is.data.frame(events_obj)) {
    # List of per-cell data.frames
    if (all(sapply(events_obj, is.data.frame))) {
      events_df <- do.call(rbind, Map(function(df, nm) {
        df$cell_id <- nm
        df
      }, events_obj, names(events_obj)))
    } else {
      stop("Cannot parse events object of class: ", class(events_obj))
    }
  } else {
    stop("Cannot parse events object of class: ", class(events_obj))
  }
  
  # Standardize column names
  names(events_df) <- tolower(names(events_df))
  # Common renames
  col_map <- c(
    "chrom" = "chr", "chromosome" = "chr", "seqnames" = "chr",
    "cell" = "cell_id", "cellid" = "cell_id", "sample" = "cell_id",
    "cn" = "copy_number", "copynumber" = "copy_number", "state" = "copy_number",
    "type" = "event_type", "event" = "event_type"
  )
  for (old_name in names(col_map)) {
    if (old_name %in% names(events_df) && !col_map[old_name] %in% names(events_df)) {
      names(events_df)[names(events_df) == old_name] <- col_map[old_name]
    }
  }
  
  # Ensure cell_id column exists
  if (!"cell_id" %in% names(events_df) && !is.null(rownames(events_df))) {
    # Check if there's a column that looks like cell IDs
    cat("  WARNING: No cell_id column found. Columns:", paste(names(events_df), collapse=", "), "\n")
  }
  
  out_events <- file.path(output_dir, "scunique_events.tsv")
  write.table(events_df, file = out_events, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  ->", nrow(events_df), "events written:", out_events, "\n")
} else {
  cat("WARNING: No uniqueEvents.RDS or df_events.RDS found\n")
}

# ─── 3. Cell metadata → CSV ──────────────────────────────────────────────────
meta_file <- find_rds("df_pass\\.RDS$")
if (is.null(meta_file)) meta_file <- find_rds("df_pass_post\\.RDS$")

if (!is.null(meta_file)) {
  cat("Loading cell metadata:", meta_file, "\n")
  meta_obj <- readRDS(meta_file)
  
  if (is.data.frame(meta_obj)) {
    meta_df <- meta_obj
    names(meta_df) <- tolower(names(meta_df))
    
    # Ensure cell_id exists
    if (!"cell_id" %in% names(meta_df)) {
      # Try common alternatives
      for (col in c("cell", "cellid", "sample", "barcode")) {
        if (col %in% names(meta_df)) {
          meta_df$cell_id <- meta_df[[col]]
          break
        }
      }
      # Last resort: use rownames
      if (!"cell_id" %in% names(meta_df)) {
        meta_df$cell_id <- rownames(meta_df)
      }
    }
    
    # Add sample_id and patient_id if missing
    if (!"sample_id" %in% names(meta_df)) meta_df$sample_id <- sample_id
    if (!"patient_id" %in% names(meta_df)) meta_df$patient_id <- sample_id
    
    out_meta <- file.path(output_dir, "cell_metadata.csv")
    write.csv(meta_df, file = out_meta, row.names = FALSE)
    cat("  ->", nrow(meta_df), "cells written:", out_meta, "\n")
  }
} else {
  cat("WARNING: No df_pass.RDS found\n")
}

# ─── Summary ─────────────────────────────────────────────────────────────────
cat("\n=== Conversion complete ===\n")
cat("Output directory:", output_dir, "\n")
cat("Files:\n")
for (f in list.files(output_dir)) cat("  ", f, "\n")

# ─── 4. BAM manifest (if bam_dir provided) ───────────────────────────────────
if (!is.null(bam_dir)) {
  cat("\nGenerating BAM manifest from:", bam_dir, "\n")
  
  # Get cell IDs from the tree
  tree_file2 <- find_rds("medicc_tree\\.RDS$")
  if (is.null(tree_file2)) tree_file2 <- find_rds("tree\\.RDS$")
  tree_obj2 <- readRDS(tree_file2)
  
  if (inherits(tree_obj2, "dendrogram")) {
    cell_ids <- labels(tree_obj2)
  } else if (inherits(tree_obj2, "phylo")) {
    cell_ids <- tree_obj2$tip.label
  } else if (inherits(tree_obj2, "hclust")) {
    cell_ids <- tree_obj2$labels
  } else {
    stop("Cannot extract cell IDs from tree of class: ", class(tree_obj2))
  }
  
  # Build manifest
  bam_manifest <- data.frame(
    cell_id    = cell_ids,
    bam_path   = file.path(bam_dir, paste0(cell_ids, ".bam")),
    bai_path   = file.path(bam_dir, paste0(cell_ids, ".bam.bai")),
    sample_id  = sample_id,
    patient_id = sample_id,
    stringsAsFactors = FALSE
  )
  
  out_manifest <- file.path(output_dir, "bam_manifest.csv")
  write.csv(bam_manifest, file = out_manifest, row.names = FALSE)
  cat("  ->", nrow(bam_manifest), "cells written:", out_manifest, "\n")
}
