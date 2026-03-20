library(data.table)

# --- Edit these as needed ---
base_dir <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep/spatial_counts_for_rctd_FIXED"
ref_counts_file <- "/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/rctd_data_prep/reference_for_rctd/sc_reference_counts.csv"

samples <- c(
  "adata_AA2__Core", "adata_Control_1", "adata_Control_2", 
  "adata_GBM1-Core", "adata_GBM1-Edge", "adata_GBM2-Core", 
  "adata_GBM2-Edge", "adata_GBM5-Edge", "adata_GBM8-Core"
)

# --- Load reference once ---
ref_counts_raw <- fread(ref_counts_file, data.table = FALSE)
rownames(ref_counts_raw) <- ref_counts_raw[, 1]
ref_counts_raw <- as.matrix(ref_counts_raw[, -1])

results <- data.frame(
  sample = samples,
  n_genes_total = NA_integer_,
  n_genes_ge3 = NA_integer_,
  frac_kept = NA_real_
)

cat("==== Checking all 9 problematic samples ====\n\n")

for (i in seq_along(samples)) {
  sample <- samples[i]
  sp_counts_file <- file.path(base_dir, paste0(sample, "_spatial_counts.csv"))
  if (!file.exists(sp_counts_file)) {
    cat(sample, ": FILE NOT FOUND\n")
    next
  }
  
  sp_counts_raw <- fread(sp_counts_file, data.table = FALSE)
  rownames(sp_counts_raw) <- sp_counts_raw[, 1]
  sp_counts_raw <- as.matrix(sp_counts_raw[, -1])
  
  # Shared genes
  shared_genes <- intersect(rownames(ref_counts_raw), rownames(sp_counts_raw))
  sp_counts <- sp_counts_raw[shared_genes, , drop = FALSE]
  
  # Count genes with >=3 counts
  gene_total_counts <- rowSums(sp_counts)
  n_genes_ge3 <- sum(gene_total_counts >= 3)
  n_genes_total <- length(gene_total_counts)
  frac_kept <- n_genes_ge3 / n_genes_total
  
  results$n_genes_total[i] <- n_genes_total
  results$n_genes_ge3[i] <- n_genes_ge3
  results$frac_kept[i] <- frac_kept
  
  cat(sample, ":\n")
  cat("   Genes with >=3 counts:", n_genes_ge3, "out of", n_genes_total, "(%.2f%%)\n", 100 * frac_kept)
  if (frac_kept < 0.1) {
    cat("   ==> WARNING: <10% of genes pass. This WILL break RCTD.\n\n")
  } else {
    cat("   OK: Sufficient genes pass threshold.\n\n")
  }
}

cat("\nSUMMARY:\n")
print(results)
