setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/new")

library(DESeq2)
library(volcano3D)
library(Seurat)
library(plotly)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(patchwork)

# ============================================================
# Load and add metadata
# ============================================================

seurat_obj <- readRDS("../merged_after_reclustering_w_res_dims20.rds")

sra1 <- read.csv("../SraRunTable (4).csv", check.names = FALSE)
sra2 <- read.csv("../SraRunTable (5).csv", check.names = FALSE)
sra  <- bind_rows(sra1, sra2)

library_df <- sra %>%
  filter(`Assay Type` == "RNA-Seq") %>%
  transmute(
    sample         = `Sample Name`,
    library_name   = `Library Name`,
    library_run    = Run,
    library_biorep = Biological_replicate,
    sra_study      = `SRA Study`
  ) %>%
  distinct()

meta  <- seurat_obj[[]]
meta2 <- meta %>%
  tibble::rownames_to_column("cell") %>%
  left_join(library_df, by = "sample") %>%
  tibble::column_to_rownames("cell")

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = meta2[, c("library_name", "library_run", "library_biorep", "sra_study")]
)

meta <- seurat_obj[[]] %>%
  mutate(
    patient_id = str_extract(library_name, "WU397-\\d{3}|WU368-\\d{2}"),
    vaccine = case_when(
      str_starts(library_name, "ELAB-WU368") ~ "covidmRNA",
      patient_id %in% c("WU397-006", "WU397-017", "WU397-028", "WU397-029") ~ "mRNA-1010",
      patient_id %in% c("WU397-005", "WU397-009", "WU397-022") ~ "Fluarix",
      TRUE ~ NA_character_
    )
  )

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = meta[, c("patient_id", "vaccine")]
)

cat("Vaccine counts:\n")
print(table(seurat_obj$vaccine, useNA = "ifany"))

# ============================================================
# Score cells by T, NK, B, Plasma and Myeloid markers
# ============================================================

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
  dplyr::group_by(RNA_snn_res.0.1) |>
  dplyr::summarise(
    t_mean    = round(mean(t_score1),     3),
    nk_mean    = round(mean(nk_score1),     3),
    b_mean      = round(mean(b_score1),       3),
    plasma_mean = round(mean(plasma_score1),  3),
    myeloid_mean = round(mean(myeloid_score1), 3),
    n_cells     = dplyr::n(),
    .groups = "drop"
  )


write.csv(score_summary, "1_score_summary.csv")

# Based on the scores,
# Define cluster-to-lineage mapping
t_clusters       <- c("0", "3", "6")
nk_clusters      <- c("4")
b_clusters       <- c("1", "7")
plasma_clusters  <- c("2")
b_plasma_clusters <- c("5")       # ambiguous — decide later
myeloid_clusters <- c("8", "9")

# Subset each lineage
t_obj       <- subset(seurat_obj, subset = RNA_snn_res.0.1 %in% t_clusters)
nk_obj      <- subset(seurat_obj, subset = RNA_snn_res.0.1 %in% nk_clusters)
b_obj       <- subset(seurat_obj, subset = RNA_snn_res.0.1 %in% b_clusters)
plasma_obj  <- subset(seurat_obj, subset = RNA_snn_res.0.1 %in% plasma_clusters)
myeloid_obj <- subset(seurat_obj, subset = RNA_snn_res.0.1 %in% myeloid_clusters)

recluster <- function(obj, resolution = 0.3, n_dims = 20, lineage_name = "lineage", out_dir = ".") {
  
  # Create output directory if it doesn't exist
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:n_dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = 1:n_dims)
  
  # ── UMAP plot ──────────────────────────────────────────────────────────────
  umap_plot <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle(paste0(lineage_name, " — res ", resolution)) +
    theme_bw()
  
  ggsave(
    filename = file.path(out_dir, paste0(lineage_name, "_umap_res", resolution, ".pdf")),
    plot     = umap_plot,
    width    = 8, height = 7
  )
  
  # ── Markers ────────────────────────────────────────────────────────────────
  all_markers <- FindAllMarkers(
    obj,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25
  )
  
  top10 <- all_markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = 10) |>
    dplyr::ungroup()
  
  write.csv(
    all_markers,
    file      = file.path(out_dir, paste0(lineage_name, "_all_markers.csv")),
    row.names = FALSE
  )
  
  write.csv(
    top10,
    file      = file.path(out_dir, paste0(lineage_name, "_top10_markers.csv")),
    row.names = FALSE
  )
  
  # ── Heatmap ────────────────────────────────────────────────────────────────
  top10_genes <- top10 |> dplyr::pull(gene) |> unique()
  
  heat_plot <- DoHeatmap(obj, features = top10_genes) +
    ggtitle(paste0(lineage_name, " — top10 markers per cluster"))
  
  ggsave(
    filename = file.path(out_dir, paste0(lineage_name, "_top10_heatmap.pdf")),
    plot     = heat_plot,
    width    = 12, height = 10
  )
  
  # ── Save RDS ───────────────────────────────────────────────────────────────
  saveRDS(
    obj,
    file = file.path(out_dir, paste0(lineage_name, "_reclustered.rds"))
  )
  
  message("[", lineage_name, "] Done — files written to: ", normalizePath(out_dir))
  
  # ── Return everything ──────────────────────────────────────────────────────
  return(list(
    obj         = obj,
    all_markers = all_markers,
    top10       = top10
  ))
}

# ── Run per lineage ───────────────────────────────────────────────────────────

out <- "results/subclustering"   # change to your preferred path

t_res       <- recluster(t_obj,       resolution = 0.4, lineage_name = "T",       out_dir = out)
nk_res      <- recluster(nk_obj,      resolution = 0.2, lineage_name = "NK",      out_dir = out)
b_res       <- recluster(b_obj,       resolution = 0.3, lineage_name = "B",       out_dir = out)
plasma_res  <- recluster(plasma_obj,  resolution = 0.2, lineage_name = "Plasma",  out_dir = out)
myeloid_res <- recluster(myeloid_obj, resolution = 0.3, lineage_name = "Myeloid", out_dir = out)


#!/usr/bin/env Rscript
# =============================================================================
# 01_convert_rds_to_h5ad.R
# Convert subsetted & reclustered Seurat RDS files to h5ad using zellkonverter
#
# Pipeline: Seurat RDS → SingleCellExperiment → writeH5AD()
# No intermediate h5Seurat file; no SeuratDisk dependency.
#
# Usage:
#   Rscript 01_convert_rds_to_h5ad.R \
#     --rds_dir  results/subclustering \
#     --out_dir  results/h5ad \
#     [--assay   RNA] \
#     [--compression gzip]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--rds_dir",
              type    = "character",
              default = "results/subclustering",
              help    = "Directory containing *_reclustered.rds files [default: %default]"),
  make_option("--out_dir",
              type    = "character",
              default = "results/h5ad",
              help    = "Output directory for .h5ad files [default: %default]"),
  make_option("--assay",
              type    = "character",
              default = "RNA",
              help    = "Seurat assay to use as AnnData X [default: %default]"),
  make_option("--compression",
              type    = "character",
              default = "gzip",
              help    = "h5ad compression: none | gzip | lzf [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Find RDS files ────────────────────────────────────────────────────────────
rds_files <- list.files(opt$rds_dir,
                        pattern    = "_reclustered\\.rds$",
                        full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No *_reclustered.rds files found in: ", opt$rds_dir)
}

message("Found ", length(rds_files), " RDS file(s):\n",
        paste(" -", basename(rds_files), collapse = "\n"))

# ── Conversion helper ─────────────────────────────────────────────────────────
convert_one <- function(rds_path, out_dir, assay = "RNA", compression = "gzip") {

  lineage <- sub("_reclustered\\.rds$", "", basename(rds_path))
  message("\n[", lineage, "] Loading RDS...")
  obj <- readRDS(rds_path)

  DefaultAssay(obj) <- assay

  # Seurat v5 stores data in layers (counts / data / scale.data).
  # as.SingleCellExperiment() handles this automatically, but we make sure
  # counts and normalised data are available as separate assay slots in SCE.
  message("[", lineage, "] Converting Seurat → SingleCellExperiment...")
  sce <- as.SingleCellExperiment(obj, assay = assay)

  # as.SingleCellExperiment maps:
  #   counts layer  → assay "counts"
  #   data layer    → assay "logcounts"   (standard SCE convention)
  #   scale.data    → assay "scaledata"   (if present)
  #   reductions    → reducedDims (PCA, UMAP, etc.)
  #   meta.data     → colData

  # Tell zellkonverter to put logcounts in X (scanpy convention after
  # normalisation); raw counts will land in a layer called "counts".
  # If you want raw counts in X instead, change X_name = "counts".
  message("[", lineage, "] Writing h5ad (X = logcounts)...")
  h5ad_path <- file.path(out_dir, paste0(lineage, ".h5ad"))

  writeH5AD(
    sce,
    file        = h5ad_path,
    X_name      = "logcounts",   # goes into AnnData.X
    compression = compression,
    verbose     = TRUE
  )

  message("[", lineage, "] Written: ", normalizePath(h5ad_path))
  return(h5ad_path)
}

# ── Run over all files ────────────────────────────────────────────────────────
h5ad_paths <- vapply(
  rds_files,
  convert_one,
  FUN.VALUE   = character(1),
  out_dir     = opt$out_dir,
  assay       = opt$assay,
  compression = opt$compression
)

message("\n=== Conversion complete ===")
message("h5ad files written to: ", normalizePath(opt$out_dir))
message(paste(" -", basename(h5ad_paths), collapse = "\n"))

