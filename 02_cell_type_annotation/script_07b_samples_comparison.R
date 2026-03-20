# Set up
library(data.table)

# ----- EDIT THIS PATH -----
counts_dir <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep/spatial_counts_for_rctd"

# List of problematic files (add or remove as needed)
problem_files <- c(
  "adata_AA2__Core_spatial_counts.csv",
  "adata_Control_1_spatial_counts.csv",
  "adata_Control_2_spatial_counts.csv",
  "adata_GBM1-Core_spatial_counts.csv",
  "adata_GBM1-Edge_spatial_counts.csv",
  "adata_GBM2-Core_spatial_counts.csv",
  "adata_GBM2-Edge_spatial_counts.csv",
  "adata_GBM5-Edge_spatial_counts.csv",
  "adata_GBM8-Core_spatial_counts.csv"
)

# Optionally, a known GOOD file for side-by-side comparison
working_file <- "adata_AA1_Core_spatial_counts.csv"

# ---- Diagnostics loop ----
cat("\n=== Checking WORKING file for reference ===\n")
ref_file <- file.path(counts_dir, working_file)
sp_work <- fread(ref_file, data.table = FALSE)
cat("\n-- WORKING FILE:", working_file, "--\n")
cat("Dimensions (rows x cols): ", dim(sp_work)[1], " x ", dim(sp_work)[2], "\n")
cat("First 5 rownames (should be gene names!):\n")
print(head(sp_work[,1], 5))
cat("First 5 column headers:\n")
print(head(colnames(sp_work), 5))
cat("\n")

for (file in problem_files) {
  cat("\n=== Problematic file:", file, "===\n")
  fpath <- file.path(counts_dir, file)
  if (!file.exists(fpath)) {
    cat("File not found! Skipping.\n")
    next
  }
  sp_prob <- fread(fpath, data.table = FALSE)
  cat("Dimensions (rows x cols): ", dim(sp_prob)[1], " x ", dim(sp_prob)[2], "\n")
  cat("First 5 rownames (should be gene names!):\n")
  print(head(sp_prob[,1], 5))
  cat("First 5 column headers:\n")
  print(head(colnames(sp_prob), 5))
  cat("First 5 data rows:\n")
  print(head(sp_prob, 5))
  cat("\n\n")
}
cat("All file diagnostics complete!\n")



#### Fixing the files#####
library(data.table)

counts_dir <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep/spatial_counts_for_rctd"
fixed_dir <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep/spatial_counts_for_rctd_FIXED"
if (!dir.exists(fixed_dir)) dir.create(fixed_dir)

problematic_samples <- c(
  "adata_AA2__Core", "adata_Control_1", "adata_Control_2",
  "adata_GBM1-Core", "adata_GBM1-Edge", "adata_GBM2-Core",
  "adata_GBM2-Edge", "adata_GBM5-Edge", "adata_GBM8-Core"
)

for (sample_name in problematic_samples) {
  cat("Processing", sample_name, "\n")
  counts_path <- file.path(counts_dir, paste0(sample_name, "_spatial_counts.csv"))
  if (!file.exists(counts_path)) {
    cat("  Missing file, skipping\n")
    next
  }
  df <- fread(counts_path, data.table = FALSE)
  # Check: is the *first column* gene names (should look like "ABCC9", "ACTA2", etc.)
  firstcol <- df[,1]
  # Heuristic: if the first value contains a dash "-", treat as barcode and force transpose
  if (grepl("-", firstcol[1])) {
    # Need to transpose
    barcodes <- df[,1]
    genes <- colnames(df)[-1]
    mat <- t(as.matrix(df[, -1]))
    colnames(mat) <- barcodes
    rownames(mat) <- genes
    df_fixed <- data.frame(Gene = rownames(mat), mat, check.names = FALSE)
    fixed_path <- file.path(fixed_dir, paste0(sample_name, "_spatial_counts.csv"))
    fwrite(df_fixed, fixed_path)
    cat("  [Fixed] File was barcodes-as-rows, now genes-as-rows. Saved to", fixed_path, "\n")
  } else {
    # Already gene names as rows
    fixed_path <- file.path(fixed_dir, paste0(sample_name, "_spatial_counts.csv"))
    fwrite(df, fixed_path)
    cat("  File seems OK (gene names in first column). Copied as is.\n")
  }
}
cat("Pre-processing done!\n")
