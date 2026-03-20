# Chiba Xenium RCTD Analysis - R Package Installation Script
# Run this script in R after setting up the Python environment

# Install CRAN packages
cran_packages <- c(
  "Matrix",
  "data.table",
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "RColorBrewer",
  "circlize",
  "devtools"
)

install.packages(cran_packages, dependencies = TRUE)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "SummarizedExperiment",
  "SpatialExperiment",
  "SingleCellExperiment",
  "S4Vectors",
  "ComplexHeatmap"
)

BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# Install spacexr (RCTD) from GitHub
# This is the core package for cell type deconvolution
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Verify installation
cat("\n=== Verifying R package installation ===\n")

all_packages <- c(cran_packages, bioc_packages, "spacexr")
installed <- sapply(all_packages, requireNamespace, quietly = TRUE)

if (all(installed)) {
  cat("\n✓ All R packages installed successfully!\n")
} else {
  cat("\n✗ Some packages failed to install:\n")
  print(all_packages[!installed])
}

# Print session info
cat("\n=== R Session Info ===\n")
sessionInfo()
