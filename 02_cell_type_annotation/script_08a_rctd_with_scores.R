# ========================== Combined RCTD runner ==========================
# A) standard pass (your original settings)
# B) fallback pass for problematic/failed samples (your tweaked settings)
# Saves: weights, scores (if found), and full .rds per sample; logs errors.
# ==========================================================================

# ---- libs ----
library(Matrix); library(data.table); library(spacexr)
library(SummarizedExperiment); library(SpatialExperiment); library(SingleCellExperiment); library(S4Vectors)

set.seed(123)

# ---- dirs ----
base_dir        <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep"
counts_dir_A    <- file.path(base_dir, "spatial_counts_for_rctd")           # Attempt A
counts_dir_B    <- file.path(base_dir, "spatial_counts_for_rctd_FIXED")     # Attempt B (if exists)
coords_dir      <- file.path(base_dir, "spatial_coords_for_rctd")
results_dir     <- file.path(base_dir, "rctd_outputs_with_scores"); if (!dir.exists(results_dir)) dir.create(results_dir)

# ---- samples known to be tricky (optional) ----
problematic_samples <- c("adata_GBM2-Edge","adata_GBM5-Edge","adata_GBM8-Core","adata_AA2__Core")

# ---- reference (shared) ----
ref_counts_file   <- file.path(base_dir, "reference_for_rctd/sc_reference_counts.csv")
ref_celltype_file <- file.path(base_dir, "reference_for_rctd/sc_reference_celltypes.csv")

ref_counts_raw <- fread(ref_counts_file, data.table = FALSE)
rownames(ref_counts_raw) <- ref_counts_raw[,1]; ref_counts_raw <- ref_counts_raw[,-1]
ref_counts_raw <- as.matrix(ref_counts_raw); storage.mode(ref_counts_raw) <- "integer"

cell_types_df <- fread(ref_celltype_file, data.table = FALSE, header = FALSE)
colnames(cell_types_df) <- c("barcode","cell_type")
cell_types_df$cell_type <- gsub("/", "_", cell_types_df$cell_type)

matched_barcodes <- intersect(colnames(ref_counts_raw), cell_types_df$barcode)
ref_counts_raw <- ref_counts_raw[, matched_barcodes, drop=FALSE]
cell_types_df  <- cell_types_df[cell_types_df$barcode %in% matched_barcodes, ]
cell_types     <- factor(cell_types_df$cell_type, levels = unique(cell_types_df$cell_type))

reference_se_full <- SummarizedExperiment(
  assays  = list(counts = ref_counts_raw),
  colData = DataFrame(cell_type = cell_types, row.names = cell_types_df$barcode)
)

# ---- helpers ----
MIN_OBS_A <- 1  # Attempt A gene_obs_min / filter
# Attempt B params (exactly your second script)
gene_cutoff_B     <- 0.000125
fc_cutoff_B       <- 0.5
gene_cutoff_reg_B <- 0.0002
fc_cutoff_reg_B   <- 0.75
UMI_min_B         <- 1
UMI_max_B         <- 20000000
UMI_min_sigma_B   <- 0
MIN_OBS_B         <- 1
max_cells_per_type_B <- 5000

# Save weights + scores + RDS
extract_and_save <- function(results_spe, sample_name) {
  weights <- as.matrix(assay(results_spe, "weights"))
  write.csv(weights, file = file.path(base_dir,   paste0("rctd_results_weights_", sample_name, ".csv")))
  write.csv(weights, file = file.path(results_dir, paste0("rctd_results_weights_", sample_name, ".csv")))

  md <- metadata(results_spe)
  scores_df <- NULL
  if ("results_df" %in% names(md)) {
    scores_df <- md$results_df
  } else if ("rctd_results" %in% names(md)) {
    scores_df <- md$rctd_results
  } else {
    cd <- as.data.frame(colData(results_spe))
    cand_names <- grep("score|likelihood|posterior|prob", names(cd), ignore.case = TRUE, value = TRUE)
    if (length(cand_names)) {
      scores_df <- cd[, cand_names, drop = FALSE]
      scores_df$pixel_id <- rownames(cd)
    }
  }
  if (!is.null(scores_df)) {
    if (is.null(scores_df$pixel_id)) scores_df$pixel_id <- rownames(scores_df)
    write.csv(scores_df, file = file.path(results_dir, paste0("rctd_results_scores_", sample_name, ".csv")), row.names = FALSE)
  } else {
    message("No explicit scores table found for ", sample_name)
  }
  saveRDS(results_spe, file = file.path(results_dir, paste0("rctd_results_object_", sample_name, ".rds")))
}

# Load counts/coords (with optional transpose check)
load_counts_coords <- function(sample_name, counts_dir) {
  sp_counts_file <- file.path(counts_dir, paste0(sample_name, "_spatial_counts.csv"))
  coords_file    <- file.path(coords_dir,  paste0(sample_name, "_spatial_coords.csv"))
  if (!file.exists(sp_counts_file)) stop("counts file missing")
  if (!file.exists(coords_file))    stop("coords file missing")

  sp <- fread(sp_counts_file, data.table = FALSE)
  rownames(sp) <- sp[,1]; sp <- sp[,-1]; sp <- as.matrix(sp)
  storage.mode(sp) <- "numeric"; sp <- round(sp); storage.mode(sp) <- "integer"

  cr <- fread(coords_file, data.table = FALSE)
  rownames(cr) <- cr[,1]; coords <- as.matrix(cr[,-1]); rownames(coords) <- cr[,1]

  if (!all(colnames(sp) %in% rownames(cr))) {
    sp <- t(sp)
    message("  Transposed spatial counts for ", sample_name)
  }
  list(counts = sp, coords = coords)
}

# Attempt A (standard)
run_attempt_A <- function(sample_name, sp_counts_raw, coords) {
  shared_genes <- intersect(rownames(assay(reference_se_full)), rownames(sp_counts_raw))
  sp_counts <- sp_counts_raw[shared_genes, , drop = FALSE]
  ref_subset <- assay(reference_se_full)[shared_genes, , drop = FALSE]

  keep <- (rowSums(sp_counts > 0) >= MIN_OBS_A) & (rowSums(ref_subset > 0) >= MIN_OBS_A)
  sp_counts  <- sp_counts[keep, , drop = FALSE]
  ref_subset <- ref_subset[keep, , drop = FALSE]

  spe <- SpatialExperiment(assays = list(counts = sp_counts), spatialCoords = coords)
  reference_se_subset <- SummarizedExperiment(assays = list(counts = ref_subset), colData = colData(reference_se_full))
  colData(reference_se_subset)$nUMI <- colSums(assay(reference_se_subset))

  rctd_obj <- createRctd(spe, reference_se_subset, ref_UMI_min = 1, UMI_min_sigma = 10, gene_obs_min = 1)
  runRctd(rctd_obj, rctd_mode = "doublet")
}

# Attempt B (tweaked) with retry ref_UMI_min 1 -> 0 and downsampling
run_attempt_B <- function(sample_name, sp_counts_raw, coords) {
  # gene filter like your second script
  genes_pass  <- rowSums(sp_counts_raw >= 1) >= MIN_OBS_B
  genes_kept  <- rownames(sp_counts_raw)[genes_pass]
  genes_final <- intersect(genes_kept, rownames(ref_counts_raw))
  if (length(genes_final) < 10) stop("Less than 10 genes left after filtering.")

  sp_counts      <- sp_counts_raw[genes_final, , drop=FALSE]
  ref_counts_sub <- ref_counts_raw[genes_final, , drop=FALSE]
  keep_cells     <- colSums(ref_counts_sub) > 0
  ref_counts_sub <- ref_counts_sub[, keep_cells, drop=FALSE]
  cell_types_sub <- cell_types_df[cell_types_df$barcode %in% colnames(ref_counts_sub), ]
  cell_types_sub <- cell_types_sub[match(colnames(ref_counts_sub), cell_types_sub$barcode), ]
  cell_types_sub$cell_type <- factor(cell_types_sub$cell_type, levels = unique(cell_types_sub$cell_type))
  if (ncol(ref_counts_sub) == 0) stop("No reference cells remain after filtering.")

  # downsample per type
  barcodes_keep <- unlist(lapply(unique(cell_types_sub$cell_type), function(ct){
    bc <- cell_types_sub$barcode[cell_types_sub$cell_type == ct]
    if (length(bc) > max_cells_per_type_B) sample(bc, max_cells_per_type_B) else bc
  }))
  ref_counts_sub <- ref_counts_sub[, barcodes_keep, drop=FALSE]
  cell_types_sub <- cell_types_sub[cell_types_sub$barcode %in% barcodes_keep, ]
  cell_types_sub <- cell_types_sub[match(colnames(ref_counts_sub), cell_types_sub$barcode), ]

  stopifnot(identical(rownames(sp_counts), rownames(ref_counts_sub)))

  spe <- SpatialExperiment(assays = list(counts = sp_counts), spatialCoords = coords)

  run_with_min <- function(ref_UMI_min_val){
    reference_se <- SummarizedExperiment(
      assays  = list(counts = ref_counts_sub),
      colData = DataFrame(cell_type = cell_types_sub$cell_type, row.names = cell_types_sub$barcode)
    )
    rctd_obj <- createRctd(
      spe, reference_se,
      gene_cutoff = gene_cutoff_B, fc_cutoff = fc_cutoff_B,
      gene_cutoff_reg = gene_cutoff_reg_B, fc_cutoff_reg = fc_cutoff_reg_B,
      ref_UMI_min = ref_UMI_min_val,
      UMI_min = UMI_min_B, UMI_max = UMI_max_B, UMI_min_sigma = UMI_min_sigma_B,
      gene_obs_min = MIN_OBS_B
    )
    runRctd(rctd_obj, rctd_mode = "doublet")
  }

  # try ref_UMI_min=1 then 0
  res <- tryCatch({
    run_with_min(1)
  }, error = function(e){
    if (grepl("no barcodes were included with nUMI at least min_UMI", conditionMessage(e))) {
      message("  Retrying with ref_UMI_min = 0 …")
      run_with_min(0)
    } else stop(e)
  })
  res
}

# ---- driver ----
all_files <- list.files(counts_dir_A, pattern = "_spatial_counts\\.csv$", full.names = TRUE)
all_samples <- gsub("_spatial_counts\\.csv$", "", basename(all_files))

successful_samples <- c()
skipped_samples    <- c()
used_attempt       <- list()
error_log          <- list()

for (sample_name in all_samples) {
  cat("\n===============================\nProcessing:", sample_name, "\n")

  # choose which counts dir to start from
  start_in_B <- sample_name %in% problematic_samples
  try_dirs <- if (start_in_B && dir.exists(counts_dir_B)) c("B","A") else c("A","B")

  result_ok <- FALSE
  for (which_dir in try_dirs) {
    counts_dir <- if (which_dir == "A") counts_dir_A else counts_dir_B
    if (!dir.exists(counts_dir)) {
      next
    }
    cat("  Attempt", which_dir, "using counts in:", counts_dir, "\n")
    # load counts/coords
    data_ok <- TRUE
    cc <- tryCatch(load_counts_coords(sample_name, counts_dir), error = function(e){ data_ok <<- FALSE; e })
    if (!data_ok) {
      cat("    Could not load data in this dir: ", conditionMessage(cc), "\n", sep = "")
      next
    }

    # run attempt
    res <- tryCatch({
      if (which_dir == "A") {
        run_attempt_A(sample_name, cc$counts, cc$coords)
      } else {
        run_attempt_B(sample_name, cc$counts, cc$coords)
      }
    }, error = function(e) e)

    if (inherits(res, "error")) {
      cat("    Attempt", which_dir, "FAILED:", conditionMessage(res), "\n")
      next
    } else {
      # success
      extract_and_save(res, sample_name)
      successful_samples <- c(successful_samples, sample_name)
      used_attempt[[sample_name]] <- which_dir
      cat("    ✔ Success with Attempt", which_dir, "\n")
      result_ok <- TRUE
      break
    }
  }

  if (!result_ok) {
    msg <- paste0("ERROR for ", sample_name, ": all attempts failed")
    skipped_samples <- c(skipped_samples, sample_name)
    error_log[[sample_name]] <- msg
    cat("  ✖ ", msg, "\n", sep = "")
  }
}

# ---- summary & log ----
cat("\nALL DONE.\n")
cat("Successful (", length(successful_samples), "): ", paste(successful_samples, collapse = ", "), "\n", sep = "")
cat("Skipped    (", length(skipped_samples), "): ", paste(skipped_samples, collapse = ", "), "\n", sep = "")

if (length(used_attempt)) {
  used <- vapply(names(used_attempt), function(s) used_attempt[[s]], character(1))
  tab  <- table(used)
  cat("Attempts used: ", paste(names(tab), tab, sep=":", collapse=", "), "\n", sep = "")
}

if (length(error_log)) {
  log_path <- file.path(results_dir, "error_log.txt")
  writeLines(unlist(error_log), con = log_path)
  cat("Wrote error log to: ", log_path, "\n", sep = "")
}
