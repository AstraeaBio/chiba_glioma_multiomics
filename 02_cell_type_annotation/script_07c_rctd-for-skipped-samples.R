# --- Load libraries ---
library(Matrix)
library(data.table)
library(spacexr)
library(SummarizedExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(S4Vectors)

set.seed(123)  # for reproducibility

# --- Set directories ---
base_dir <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep"
counts_dir <- file.path(base_dir, "spatial_counts_for_rctd_FIXED")
coords_dir <- file.path(base_dir, "spatial_coords_for_rctd")
results_dir <- file.path(base_dir, "rctd_outputs")
if (!dir.exists(results_dir)) dir.create(results_dir)

# --- Sample list ---
problematic_samples <- c(
  "adata_GBM2-Edge",
  "adata_GBM5-Edge", "adata_GBM8-Core", "adata_AA2__Core",
)

# --- Reference files ---
ref_counts_file <- file.path(base_dir, "reference_for_rctd/sc_reference_counts.csv") # nolint # nolint
ref_celltype_file <- file.path(base_dir, "reference_for_rctd/sc_reference_celltypes.csv") # nolint: line_length_linter.

# --- Load reference (raw, unfiltered) ---
ref_counts_raw <- fread(ref_counts_file, data.table = FALSE)
rownames(ref_counts_raw) <- ref_counts_raw[, 1]
ref_counts_raw <- as.matrix(ref_counts_raw[, -1])
storage.mode(ref_counts_raw) <- "integer"

cell_types_df <- fread(ref_celltype_file, data.table = FALSE, header = FALSE)
colnames(cell_types_df) <- c("barcode", "cell_type")
cell_types_df$cell_type <- gsub("/", "_", cell_types_df$cell_type)
matched_barcodes <- intersect(colnames(ref_counts_raw), cell_types_df$barcode)
ref_counts_raw <- ref_counts_raw[, matched_barcodes, drop = FALSE]
cell_types_df <- cell_types_df[cell_types_df$barcode %in% matched_barcodes, ]

# --- RCTD config parameters ---
gene_cutoff <- 0.000125
fc_cutoff <- 0.5
gene_cutoff_reg <- 0.0002
fc_cutoff_reg <- 0.75
UMI_min <- 1 # nolint
UMI_max <- 20000000 # nolint 
UMI_min_sigma <- 0 # nolint
MIN_OBS <- 1    # for sample filtering and gene_obs_min in createRctd # nolint

# --- Loop over samples ---
for (sample in problematic_samples) {
  cat("\n===============================\nProcessing:", sample, "\n")
   # nolint
  sp_counts_file <- file.path(counts_dir, paste0(sample, "_spatial_counts.csv"))
  coords_file    <- file.path(coords_dir, paste0(sample, "_spatial_coords.csv"))
  if (!file.exists(sp_counts_file) | !file.exists(coords_file)) { # nolint
    cat("  Missing files for", sample, "; skipping.\n")
    next
  }

  # --- Load spatial counts ---
  sp_counts_raw <- fread(sp_counts_file, data.table = FALSE)
  rownames(sp_counts_raw) <- sp_counts_raw[, 1]
  sp_counts_raw <- as.matrix(sp_counts_raw[, -1])
  storage.mode(sp_counts_raw) <- "integer"
   # nolint # nolint
  # --- Load coordinates ---
  coords_raw <- fread(coords_file, data.table = FALSE)
  rownames(coords_raw) <- coords_raw[, 1]
  coords <- as.matrix(coords_raw[, -1])
  rownames(coords) <- coords_raw[, 1]
   # nolint # nolint
  # --- Per-sample gene filter ---
  genes_pass <- rowSums(sp_counts_raw >= 1) >= MIN_OBS
  genes_kept <- rownames(sp_counts_raw)[genes_pass]
  genes_ref  <- rownames(ref_counts_raw)
  genes_final <- intersect(genes_kept, genes_ref)
  if (length(genes_final) < 10) {
    cat("  WARNING: Less than 10 genes left after filtering. Skipping sample.\n") # nolint
    next
  }
   # nolint
  sp_counts <- sp_counts_raw[genes_final, , drop = FALSE]
  ref_counts_sub <- ref_counts_raw[genes_final, , drop = FALSE]
  keep_cells <- colSums(ref_counts_sub) > 0
  ref_counts_sub <- ref_counts_sub[, keep_cells, drop = FALSE]
  cell_types_sub <- cell_types_df[cell_types_df$barcode %in% colnames(ref_counts_sub), ] # nolint
  cell_types_sub <- cell_types_sub[match(colnames(ref_counts_sub), cell_types_sub$barcode), ] # nolint # nolint
  cell_types_sub$cell_type <- factor(cell_types_sub$cell_type, levels = unique(cell_types_sub$cell_type)) # nolint

  n_ref_cells <- ncol(ref_counts_sub)
  min_UMI_obs <- ifelse(n_ref_cells > 0, min(colSums(ref_counts_sub)), NA) # nolint
  cat("  Reference cells after gene/cell filtering:", n_ref_cells, "\n")
  cat("  Min UMI in reference cells:", min_UMI_obs, "\n")
   # nolint
  if (n_ref_cells == 0) {
    cat("  No reference cells remain after filtering, skipping sample.\n")
    next
  }
   # nolint
  # --- Create SpatialExperiment & reference SummarizedExperiment ---
  spe <- SpatialExperiment(
    assays = list(counts = sp_counts),
    spatialCoords = coords
  )
  stopifnot(identical(colnames(ref_counts_sub), cell_types_sub$barcode))
  stopifnot(identical(rownames(ref_counts_sub), rownames(spe)))

   # nolint
cat("Check dim/sp_counts:", dim(sp_counts), "\n") # nolint
cat("Check dim/ref_counts_sub:", dim(ref_counts_sub), "\n")
cat("All genes identical? ", identical(rownames(sp_counts), rownames(ref_counts_sub)), "\n") # nolint # nolint
cat("All barcodes identical? ", identical(colnames(ref_counts_sub), cell_types_sub$barcode), "\n") # nolint


# Double-check # nolint
stopifnot(identical(rownames(sp_counts), rownames(ref_counts_sub)))
stopifnot(identical(colnames(ref_counts_sub), cell_types_sub$barcode))
print(table(cell_types_sub$cell_type))
print(rownames(sp_counts)[!rownames(sp_counts) %in% rownames(ref_counts_sub)])


# --- Manual downsampling per cell type --- # nolint
max_cells_per_type <- 5000
barcodes_keep <- c()
for (ct in unique(cell_types_sub$cell_type)) {
  these_barcodes <- cell_types_sub$barcode[cell_types_sub$cell_type == ct]
  if (length(these_barcodes) > max_cells_per_type) {
    barcodes_keep <- c(barcodes_keep, sample(these_barcodes, max_cells_per_type)) # nolint
  } else {
    barcodes_keep <- c(barcodes_keep, these_barcodes)
  }
}
ref_counts_sub <- ref_counts_sub[, barcodes_keep, drop=FALSE] # nolint # nolint
cell_types_sub <- cell_types_sub[cell_types_sub$barcode %in% barcodes_keep, ]
cell_types_sub <- cell_types_sub[match(colnames(ref_counts_sub), cell_types_sub$barcode), ] # nolint
cat("Reference after manual downsampling:\n")

print(table(cell_types_sub$cell_type)) # nolint # nolint

cat("New reference dim:", dim(ref_counts_sub), "\n") # nolint
cat("sp_counts dim: ", dim(sp_counts), "\n")
cat("ref_counts_sub dim: ", dim(ref_counts_sub), "\n")
cat("All genes identical? ", identical(rownames(sp_counts), rownames(ref_counts_sub)), "\n") # nolint
cat("All reference barcodes identical to cell_types_sub? ", identical(colnames(ref_counts_sub), cell_types_sub$barcode), "\n") # nolint
cat("Are all reference barcodes unique? ", length(unique(colnames(ref_counts_sub))) == ncol(ref_counts_sub), "\n") # nolint
cat("Number of reference cell types: ", length(unique(cell_types_sub$cell_type)), "\n") # nolint


ct_table <- table(cell_types_sub$cell_type) # nolint # nolint
cat("Cell type counts just before RCTD:\n")
print(ct_table)
if(any(ct_table == 0)) stop("At least one cell type has zero cells!") # nolint

reference_se <- SummarizedExperiment( # nolint
    assays = list(counts = ref_counts_sub),
    colData = DataFrame(cell_type = cell_types_sub$cell_type, row.names = cell_types_sub$barcode) # nolint
  )

  # --- RCTD with full config arguments ---
  run_rctd <- function(ref_UMI_min_val) { # nolint
    createRctd(
      spe, reference_se, # nolint
      gene_cutoff = gene_cutoff,
      fc_cutoff = fc_cutoff,
      gene_cutoff_reg = gene_cutoff_reg,
      fc_cutoff_reg = fc_cutoff_reg,
      ref_UMI_min = ref_UMI_min_val,
      UMI_min = UMI_min,
      UMI_max = UMI_max,
      UMI_min_sigma = UMI_min_sigma,
      gene_obs_min = MIN_OBS)
  }
   # nolint
  # --- Try running RCTD, retry with ref_UMI_min = 0 if needed ---
  rctd_status <- tryCatch({
    rctd_data <- run_rctd(1)
    cat("    RCTD object created. Running deconvolution...\n")
    results_spe <- runRctd(rctd_data, rctd_mode = "doublet")
    weights <- as.matrix(assay(results_spe, "weights"))
    out_file <- file.path(results_dir, paste0("rctd_results_weights_", sample, ".csv")) # nolint
    write.csv(weights, file = out_file)
    cat("    Results saved to", out_file, "\n")
    "success"
  }, error = function(e) {
    msg <- conditionMessage(e)
    cat("    ERROR for", sample, ":", msg, "\n")
    if (grepl("no barcodes were included with nUMI at least min_UMI", msg)) {
      cat("    Lowering ref_UMI_min to 0 and retrying...\n")
      tryCatch({
        rctd_data <- run_rctd(0)
        results_spe <- runRctd(rctd_data, rctd_mode = "doublet")
        weights <- as.matrix(assay(results_spe, "weights"))
        out_file <- file.path(results_dir, paste0("rctd_results_weights_", sample, ".csv")) # nolint
        write.csv(weights, file = out_file)
        cat("    Results saved to", out_file, "\n")
        "success_retry"
      }, error = function(e2) {
        cat("    ERROR after retry for", sample, ":", conditionMessage(e2), "\n") # nolint
        "failed"
      })
    } else {
      "failed"
    }
  })
   # nolint # nolint
  cat("    Status for", sample, ":", rctd_status, "\n")
}
