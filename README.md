# Integrated Single Cell Spatial Multi-Omics Landscape of WHO Grades 2–4 Gliomas

## Code Repository

This repository contains the analysis code for the manuscript:

> Ma Y, Ayyadhury S, Singh S, Vashishath Y, Ozdemir C, Liebenberg K, McKee TD, Nguyen N, Basi A, Mak D, Gomez JA, Huse JT, Noor S, Winkowski D, Baird R, Lang FF, Burks JK, Bozdag S, Weinberg JS, Seeley EH, Eberlin LS, Ene CI. "Integrated single cell spatial multi-omics landscape of WHO grades 2–4 gliomas identifies locoregional metabolomic surrogates of glioma transcriptional cellular states." (2026).

---

## Repository Structure

Scripts are organized by analysis stage and mapped to manuscript figures.

| Folder | Description | Figures |
|--------|-------------|---------|
| `01_preprocessing/` | Xenium QC, image registration (VALIS) | Supp. Fig. 1 |
| `02_cell_type_annotation/` | GBmap label transfer, RCTD deconvolution | Supp. Fig. 2 |
| `03_figure2_cellular_heterogeneity/` | Cell type proportions, heatmaps, UMAPs, statistical tests | Fig. 2, Supp. Fig. 3 |
| `04_figure3_deg_survival/` | Differential gene expression (Moran's I, Wilcoxon), survival analysis | Fig. 3 |
| `05_figure4_spatial_interactions/` | SCIMAP spatial cell-cell interaction analysis | Fig. 4 |
| `06_figure5_tgfb_pathway/` | TGFβ pathway expression analysis and statistical comparisons | Fig. 5 |
| `07_figure6_msi_integration/` | MSI-Xenium co-registration, metabolite-cell state integration, Ki-67 analysis | Fig. 6 |
| `08_imc_analysis/` | Imaging mass cytometry segmentation, quantification, registration | Supp. Figs. 4–5 |
| `config/` | Environment files, R package installation | — |

---

## Figure-to-Script Map

### Figure 1 — Workflow Schematic
No analysis code. Created with BioRender.com.

### Figure 2 — Cellular Heterogeneity
| Panel | Script | Language |
|-------|--------|----------|
| 2a (PCA) | `02_cell_type_annotation/script_08b_rctd_label-transfer.ipynb` | Python |
| 2b (stacked bar) | `03_figure2_cellular_heterogeneity/script_04-analysis-celltype-proportions.ipynb` | Python |
| 2c–d (heatmaps) | `03_figure2_cellular_heterogeneity/script_08c_heatmaps-for-figure2.R` | R |
| 2e–g (UMAPs) | `03_figure2_cellular_heterogeneity/script_08d_revising_figures2.ipynb` | Python |
| 2h–j (box plots + stats) | `03_figure2_cellular_heterogeneity/proportion_statistical_tests.py` | Python |

### Figure 3 — DEG and Survival
| Panel | Script | Language |
|-------|--------|----------|
| 3a–c (DEG heatmaps) | `04_figure3_deg_survival/script_05-analysis_top_ranked_analysis_merged_groups.ipynb` | Python |
| 3a–c (Moran's I) | `04_figure3_deg_survival/script_09a_rctd_spatially_variable_genes_moransI.ipynb` | Python |
| 3d–g (survival) | `04_figure3_deg_survival/survival_analysis.R` | R |

### Figure 4 — Cell-Cell Interactions
| Panel | Script | Language |
|-------|--------|----------|
| 4a–b (Oligo) | `05_figure4_spatial_interactions/script_09b_rctd_neighbourhood_analysis_ol.ipynb` | Python |
| 4c–d (AA) | `05_figure4_spatial_interactions/script_09b_rctd_neighbourhood_analysis_aa.ipynb` | Python |
| 4e–f (GBM) | `05_figure4_spatial_interactions/script_09b_rctd_neighbourhood_analysis_gbm.ipynb` | Python |

### Figure 5 — TGFβ Pathway
| Panel | Script | Language |
|-------|--------|----------|
| 5a–e (all) | `06_figure5_tgfb_pathway/Chiba_TGFB_Figure_Tumor_vs_Immune.ipynb` | Python |

### Figure 6 — MSI Integration
| Panel | Script | Language |
|-------|--------|----------|
| 6b–c (Spearman heatmaps) | `07_figure6_msi_integration/script01_analysis.ipynb` | Python |
| 6d (Ki-67 enrichment) | `07_figure6_msi_integration/P3_07_mki67_ido1_analysis.ipynb` | Python |
| 6e (MSPen) | `07_figure6_msi_integration/mspen_validation.ipynb` | Python |
| MSI co-registration | `07_figure6_msi_integration/add_warped_coords_to_msi_adatas_v4_CORRECT.py` | Python |

### Supplementary Figures
| Figure | Script |
|--------|--------|
| Supp. Fig. 1 (QC) | `01_preprocessing/script01c_intitial_processing_adatas.ipynb` |
| Supp. Fig. 2 (GBmap) | `02_cell_type_annotation/script_03-gbmap_annotation.ipynb` |
| Supp. Fig. 3 (proportions) | `03_figure2_cellular_heterogeneity/proportion_statistical_tests.py` |
| Supp. Fig. 4–5 (IMC) | `08_imc_analysis/script01_segment_IMC_cores_extract_protein_expression.ipynb` |

---

## System Requirements

### Python (≥ 3.10)
Core packages: scanpy, anndata, squidpy, scimap (≥ 2.3.5), harmonypy, scipy, statsmodels, pandas, numpy, matplotlib, seaborn

### R (≥ 4.3)
Core packages: spacexr (v1.0.0), survival, survminer, ComplexHeatmap, SpatialExperiment, circlize

### Environment Setup

```bash
# Create conda environment
conda env create -f config/environment.yml
conda activate chiba_xenium_rctd

# Install R packages (run from within the conda environment)
Rscript config/install_r_packages.R
```

---

## Data Availability

Raw and processed data are deposited in the National Metabolomics Data Repository and will be made publicly available upon publication. The tissue microarray layout and clinical information are provided in Supplementary Tables 1–2 of the manuscript.

---

## License

This code is made available for academic use under the MIT License. Please cite the manuscript if you use any portion of this code.

## Contact

Chibawanye I. Ene — cene@mdanderson.org  
The University of Texas MD Anderson Cancer Center, Houston TX 77030

For technical questions, contact
Trevor D. McKee - tdmckee@alum.mit.edu
