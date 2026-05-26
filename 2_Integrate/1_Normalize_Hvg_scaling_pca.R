#!/usr/bin/env Rscript

setwd("/media/csbl/sandbox-SSD-3/mra_vaccine")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

dir.create("results", showWarnings = FALSE)
dir.create("results/metadata_audit", showWarnings = FALSE)

# =============================================================
# Load merged raw object
# =============================================================

message("Loading merged raw object...")
merged <- readRDS("results/merged_raw.rds")

message("  Cells: ", ncol(merged))
message("  Genes: ", nrow(merged))

# =============================================================
# RNA assay
# =============================================================

DefaultAssay(merged) <- "RNA"

# save
saveRDS(
  merged,
  "results/merged_loaded.rds"
)

# =============================================================
# Normalize
# =============================================================

message("\nNormalizing...")

merged <- NormalizeData(
  merged,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

saveRDS(
  merged,
  "results/merged_normalized.rds"
)

message("  Saved: merged_normalized.rds")

# =============================================================
# HVG selection
# =============================================================

message("\nFinding HVGs...")

merged <- FindVariableFeatures(
  merged,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

message(
  "  HVGs selected: ",
  length(VariableFeatures(merged))
)

saveRDS(
  merged,
  "results/merged_hvg.rds"
)

message("  Saved: merged_hvg.rds")

# =============================================================
# Scaling
# =============================================================

message("\nScaling data...")

merged <- ScaleData(
  merged,
  features = VariableFeatures(merged),
  vars.to.regress = "percent_mt",
  verbose = FALSE
)

saveRDS(
  merged,
  "results/merged_scaled.rds"
)

message("  Saved: merged_scaled.rds")

# =============================================================
# PCA
# =============================================================

message("\nRunning PCA...")

merged <- RunPCA(
  merged,
  features = VariableFeatures(merged),
  npcs = 50,
  verbose = FALSE
)

saveRDS(
  merged,
  "results/merged_pca.rds"
)

message("  Saved: merged_pca.rds")

# =============================================================
# Elbow plot
# =============================================================

png(
  "results/metadata_audit/elbow_plot.png",
  width = 800,
  height = 500
)

ElbowPlot(
  merged,
  ndims = 50
)

dev.off()

message("  Saved elbow plot")

message("\nDone.")
