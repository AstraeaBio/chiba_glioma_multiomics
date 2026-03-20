library(Matrix)
library(data.table)
library(spacexr)
library(SummarizedExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(S4Vectors)

base_dir <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep"
counts_dir <- file.path(base_dir, "spatial_counts_for_rctd")
coords_dir <- file.path(base_dir, "spatial_coords_for_rctd")

# ---------- Define/Make results directory ----------
results_dir <- file.path(base_dir, "rctd_outputs")
if (!dir.exists(results_dir)) dir.create(results_dir)

# ----------- REFERENCE: LOAD AS-IS ------------
ref_counts_file <- file.path(base_dir, "reference_for_rctd/sc_reference_counts.csv")
ref_celltype_file <- file.path(base_dir, "reference_for_rctd/sc_reference_celltypes.csv")

cat("Loading reference counts...\n")
ref_counts_raw <- fread(ref_counts_file, data.table = FALSE)
rownames(ref_counts_raw) <- ref_counts_raw[, 1]
ref_counts_raw <- ref_counts_raw[, -1]
ref_counts_raw <- as.matrix(ref_counts_raw)
storage.mode(ref_counts_raw) <- "integer"

cat("Loading cell type labels...\n")
cell_types_df <- fread(ref_celltype_file, data.table = FALSE, header = FALSE)
colnames(cell_types_df) <- c("barcode", "cell_type")
cat("Counts shape (raw):", dim(ref_counts_raw)[1], "genes x", dim(ref_counts_raw)[2], "cells\n")

matched_barcodes <- intersect(colnames(ref_counts_raw), cell_types_df$barcode)
cat("Matched barcodes:", length(matched_barcodes), "\n")
ref_counts <- ref_counts_raw[, matched_barcodes, drop = FALSE]
cell_types_df <- cell_types_df[cell_types_df$barcode %in% matched_barcodes, ]
cell_types_df$cell_type <- gsub("/", "_", cell_types_df$cell_type)
cell_types <- factor(cell_types_df$cell_type, levels = unique(cell_types_df$cell_type))
cat("Levels in cell_types factor (after cleaning):\n")
print(levels(cell_types))

reference_se <- SummarizedExperiment(
    assays = list(counts = ref_counts),
    colData = DataFrame(cell_type = cell_types, row.names = cell_types_df$barcode)
)
cat("Reference SummarizedExperiment created: ", dim(assay(reference_se)), "\n")
umi_per_cell <- colSums(assay(reference_se))
cat("Summary of UMI per cell:\n")
print(summary(umi_per_cell))

# -------- LOOP OVER ALL SAMPLES IN counts_dir -----------
spatial_files <- list.files(counts_dir, pattern = "_spatial_counts\\.csv$", full.names = TRUE)
cat("Found", length(spatial_files), "spatial counts files.\n")
MIN_OBS <- 1  # You can set to 2 or more if needed

skipped_samples <- c()
successful_samples <- c()

for (sp_counts_file in spatial_files) {
    sample_name <- gsub("_spatial_counts\\.csv$", "", basename(sp_counts_file))
    cat("\nProcessing sample:", sample_name, "\n")
    error_flag <- FALSE

    # Load spatial counts
    sp_counts_raw <- fread(sp_counts_file, data.table = FALSE)
    rownames(sp_counts_raw) <- sp_counts_raw[, 1]
    sp_counts_raw <- sp_counts_raw[, -1]
    sp_counts_raw <- as.matrix(sp_counts_raw)
    storage.mode(sp_counts_raw) <- "numeric"
    sp_counts_raw <- round(sp_counts_raw)
    storage.mode(sp_counts_raw) <- "integer"

    # Load coords
    coords_file <- file.path(coords_dir, paste0(sample_name, "_spatial_coords.csv"))
    if (!file.exists(coords_file)) {
        cat("  No coords file for", sample_name, "; skipping.\n")
        next
    }
    coords_raw <- fread(coords_file, data.table = FALSE)
    rownames(coords_raw) <- coords_raw[, 1]
    coords <- as.matrix(coords_raw[, -1])
    rownames(coords) <- coords_raw[, 1]

    # Transpose spatial counts if needed
    if (!all(colnames(sp_counts_raw) %in% rownames(coords_raw))) {
        sp_counts_raw <- t(sp_counts_raw)
        cat("  Transposed spatial counts for", sample_name, "\n")
    }

    # Calculate shared genes
    shared_genes <- intersect(rownames(ref_counts), rownames(sp_counts_raw))
    cat("Reference genes:", head(rownames(ref_counts)), "\n")
    cat("Spatial genes:", head(rownames(sp_counts_raw)), "\n")
    cat("Number of shared genes:", length(shared_genes), "\n")
    cat("Example shared genes:", head(shared_genes), "\n")

    # Subset to shared genes
    sp_counts <- sp_counts_raw[shared_genes, , drop = FALSE]
    ref_counts_subset <- ref_counts[shared_genes, , drop = FALSE]

    # Filter out genes all-zero in either matrix
    nonzero_in_spatial <- rowSums(sp_counts > 0) >= MIN_OBS
    nonzero_in_reference <- rowSums(ref_counts_subset > 0) >= MIN_OBS
    keep_genes <- nonzero_in_spatial & nonzero_in_reference

    cat("  Genes with >=", MIN_OBS, "nonzero in spatial: ", sum(nonzero_in_spatial), "\n")
    cat("  Genes with >=", MIN_OBS, "nonzero in reference: ", sum(nonzero_in_reference), "\n")
    cat("  Genes passing both: ", sum(keep_genes), "\n")
    cat("  (First 10 kept genes): ", paste(head(rownames(sp_counts)[keep_genes]), collapse = ", "), "\n")

    if(sum(keep_genes) < 20) {
        cat("  WARNING: Very few genes passing filter; RCTD may fail for this sample.\n")
    }

    # Final filtered counts
    sp_counts <- sp_counts[keep_genes, , drop = FALSE]
    ref_counts_final <- ref_counts_subset[keep_genes, , drop = FALSE]

    # Build/Update SpatialExperiment
    spe <- SpatialExperiment(
        assays = list(counts = sp_counts),
        spatialCoords = coords
    )
    cat("  SpatialExperiment object created for", sample_name, "\n")
    cat("    genes x pixels after filtering:", dim(sp_counts), "\n")

    # Build a reference SummarizedExperiment for this sample
    reference_se_subset <- SummarizedExperiment(
        assays = list(counts = ref_counts_final),
        colData = colData(reference_se)
    )
    ref_nUMI <- colSums(assay(reference_se_subset))
    colData(reference_se_subset)$nUMI <- ref_nUMI

    # Diagnostics
    cat("Summary of nUMI per cell:\n")
    print(summary(ref_nUMI))
    cat("Number of cells with UMI > 0: ", sum(ref_nUMI > 0), "\n")
    zero_cells <- which(ref_nUMI == 0)
    cat("Cells with all-zero UMI: ", length(zero_cells), "\n")
    if (length(zero_cells) > 0) {
        cat("First 5 all-zero cell barcodes: ", head(colnames(assay(reference_se_subset))[zero_cells]), "\n")
    }
    cat("Sample:", sample_name, "\n")
    cat("Num shared genes:", length(shared_genes), "\n")
    cat("sp_counts all zeros per gene:", sum(rowSums(sp_counts) == 0), "\n")
    cat("ref_counts all zeros per gene:", sum(rowSums(ref_counts_final) == 0), "\n")

    # --- RCTD WITH ERROR HANDLING ---
    result <- tryCatch({
        rctd_data <- createRctd(spe, reference_se_subset, ref_UMI_min = 1, UMI_min_sigma = 10, gene_obs_min = 1)
        cat("  RCTD object created. Running deconvolution...\n")
        results_spe <- runRctd(rctd_data, rctd_mode = "doublet")
        cat("  runRctd complete for", sample_name, "\n")
        # Save results to *both* main and outputs folder
        weights <- as.matrix(assay(results_spe, "weights"))
        out_file_main <- file.path(base_dir, paste0("rctd_results_weights_", sample_name, ".csv"))
        out_file_sep  <- file.path(results_dir, paste0("rctd_results_weights_", sample_name, ".csv"))
        write.csv(weights, file = out_file_main)
        write.csv(weights, file = out_file_sep)
        cat("  Results saved to", out_file_main, "and", out_file_sep, "\n")
        successful_samples <<- c(successful_samples, sample_name)
        TRUE
    }, error = function(e) {
        cat("  ERROR for", sample_name, ":", conditionMessage(e), "\n")
        skipped_samples <<- c(skipped_samples, sample_name)
        FALSE
    })
}

cat("\nALL DONE!\n")
cat("Files with errors/skipped:\n")
print(skipped_samples)
cat("Files completed successfully:\n")
print(successful_samples)
