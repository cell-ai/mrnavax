setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

# ─────────────────────────────────────────────
# Load object and add cell-type annotation
# ─────────────────────────────────────────────

seurat_obj <- readRDS("merged_after_reclustering_no_TCR_BCR.rds")

t_markers <- list(c("CD3E", "CD4", "CD8A", "CD8B",
                    "FOXP3", "IL2RA",               # Treg
                    "BCL6", "CXCR5",                # Tfh
                    "CCR7", "SELL",                 # naive/central memory
                    "GZMK", "GZMB"))                # effector/cytotoxic T

nk_markers <- list(c("GNLY", "NKG7", "NCAM1",      # pan-NK
                     "KLRF1", "KLRD1",              # NK receptors
                     "FCGR3A",                      # CD16 — mature NK
                     "XCL1", "XCL2"))               # NK regulatory subset

b_markers    <- list(c("CD19", "MS4A1", "CD79A", "CD27", "IGHD", "AICDA",
                       "BCL6", "TCL1A"))

myeloid_markers <- list(
  c("CD14", "LYZ", "S100A8",        # Classical monocytes
    "FCGR3A", "MS4A7",              # Non-classical monocytes
    "LILRA4", "IRF7", "CLEC4C",     # pDC
    "FCER1A", "CD1C",               # cDC2
    "CLEC9A", "XCR1",               # cDC1
    "CD68", "MRC1")                 # Macrophages
)


plasma_markers <- list(c("MZB1", "XBP1", "PRDM1", "CD38", "IGHA1", "IGHG1"))

seurat_obj <- AddModuleScore(seurat_obj, features = t_markers,  name = "t_score")
seurat_obj <- AddModuleScore(seurat_obj, features = nk_markers, name = "nk_score")
seurat_obj <- AddModuleScore(seurat_obj, features = b_markers,       name = "b_score")
seurat_obj <- AddModuleScore(seurat_obj, features = plasma_markers,  name = "plasma_score")
seurat_obj <- AddModuleScore(seurat_obj, features = myeloid_markers, name = "myeloid_score")

# Summary table with all four scores
score_summary <- seurat_obj@meta.data |>
  dplyr::group_by(RNA_snn_res.0.1, celltype_l1) |>
  dplyr::summarise(
    t_mean    = round(mean(t_score1),     3),
    nk_mean    = round(mean(nk_score1),     3),
    b_mean      = round(mean(b_score1),       3),
    plasma_mean = round(mean(plasma_score1),  3),
    myeloid_mean = round(mean(myeloid_score1), 3),
    n_cells     = dplyr::n(),
    .groups = "drop"
  )


print(score_summary)

# UMAP with all four scores
FeaturePlot(seurat_obj,
            features  = c("t_score1", "nk_score1", "b_score1",
                          "plasma_score1", "myeloid_score1"),
            ncol      = 2,
            order     = TRUE,
            pt.size   = 0.3) &
  scale_colour_viridis_c(option = "magma") &
  theme_minimal(base_size = 10)

saveRDS(seurat_obj, "atual_seurat_obj_notcrbcr_forrecluster.rds")
