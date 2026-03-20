#install.packages(c("circlize","readr","dplyr"), repos="https://cloud.r-project.org")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(readr)
library(dplyr)

# ---------- helpers ----------
read_matrix <- function(mat_csv) {
  df <- readr::read_csv(mat_csv, show_col_types = FALSE)
  rn <- df[[1]]
  mat <- as.matrix(df[,-1,drop=FALSE])
  rownames(mat) <- rn
  mat
}
read_rowsplit <- function(csv) {
  tib <- readr::read_csv(csv, show_col_types = FALSE)
  rs <- tib$block; names(rs) <- tib$gene; rs
}
hm_colors <- circlize::colorRamp2(c(-2,0,2), c("#0ca3ae","#f7f7f7","#b218af"))  # similar to 'vlag'

# group colors
macro_cols <- c("Tumor Type"="#3a5568","BV"="#9467bd","Immune"="#e1573b")

# ---------- LEFT panel (3 blocks only) ----------
Z_left <- read_matrix("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/complexheatmap_inputs/Z_left.csv")        
row_split_left <- read_rowsplit("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/complexheatmap_inputs/Z_left_rowSplit.csv")

# Ensure ordering
row_split_left <- factor(row_split_left[rownames(Z_left)],
                         levels = c("Tumor Type","BV","Immune"))

# --- annotations ---
col_annot <- HeatmapAnnotation(
  Group = factor(colnames(Z_left), levels=c("Tumor Type","BV","Immune")),
  col = list(Group = c(
    "Tumor Type"="#3a5568",
    "BV"="#9467bd","Immune"="#e1573b")),
  annotation_name_side = "left"
)

row_annot <- rowAnnotation(
  Block = row_split_left,
  col = list(Block = macro_cols),
  annotation_name_side = "top", width = unit(6, "mm")
)

# --- heatmap ---
ht_left <- Heatmap(
  Z_left,
  name = "z-score",
  col = hm_colors,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_split = row_split_left,     # ✅ only 3 blocks now
  top_annotation = col_annot,
  left_annotation = row_annot,
  column_title = "Tumor / BV / Immune - top markers",
  heatmap_legend_param = list(title = "z-score"),
  border = TRUE
)

# save
pdf("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/fig2c_left_complexheatmap.pdf", width=7.2, height=6)
draw(ht_left, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/fig2c_left_complexheatmap.png", width=2500, height=2500, res=300)
draw(ht_left, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# ---------- Immune-only panel (refined) ----------
# choose how many markers per subtype
top_n <- 3   # set to 3 or 9 as you prefer

# read immune subtype matrix
Z_imm <- read_matrix("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/complexheatmap_inputs/Z_imm.csv")
rsI <- read_rowsplit("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/complexheatmap_inputs/Z_imm_rowSplit.csv")
rsI <- factor(rsI[rownames(Z_imm)], levels=unique(rsI[rownames(Z_imm)]))

# ---- filter top N genes per block ----
keep_genes <- unlist(
  tapply(rownames(Z_imm), rsI[rownames(Z_imm)], function(x) head(x, top_n))
)
Z_imm_sub <- Z_imm[keep_genes, ]
rsI_sub <- rsI[keep_genes]

# ---- subtype colors ----
subtype_names <- colnames(Z_imm_sub)
subtype_palette <- RColorBrewer::brewer.pal(n = length(subtype_names), name = "Set3")
names(subtype_palette) <- subtype_names

col_annot_I <- HeatmapAnnotation(
  Subtype = factor(colnames(Z_imm_sub), levels=subtype_names),
  col = list(Subtype = subtype_palette),
  annotation_name_side = "left"
)

row_annot_I <- rowAnnotation(
  Block = rsI_sub,
  col = list(Block = subtype_palette[levels(rsI_sub)]),
  annotation_name_side = "top", width = unit(8, "mm")
)

ht_imm <- Heatmap(
  Z_imm_sub,
  name = "z-score",
  col = hm_colors,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_split = rsI_sub,
  top_annotation = col_annot_I,
  left_annotation = row_annot_I,
  column_title = paste("Immune subtypes — top", top_n, "markers"),
  column_names_gp = gpar(fontsize=10, fontface="bold"),
  row_names_gp = gpar(fontsize=8),
  row_title_gp = gpar(fontsize=10, fontface="bold"),
  row_title_rot = 0,
  border = TRUE,
  heatmap_legend_param = list(
    title = "z-score",
    legend_height = unit(3, "cm"),
    labels_gp = gpar(fontsize=8)
  )
)

pdf("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/fig2c_right_complexheatmap.pdf",
    width = 8.2, height = 6)
draw(ht_imm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

png("/home/ext_sana_noor_astraeabio_com/ext_hd/chiba/adatas_with_rctd/figs/fig2c_right_complexheatmap.png",
    width = 2300, height = 1500, res = 300)
draw(ht_imm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()


