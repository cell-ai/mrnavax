# =============================================================
# Harmony Integration
# =============================================================
message("\nRunning Harmony integration...")

library(harmony)

merged <- readRDS("results/merged_pca.rds")

# =============================================================
# Check / recreate 'study' metadata column from cell ID prefix
# =============================================================
if (!"study" %in% colnames(merged@meta.data)) {
  message("  'study' column not found — recreating from cell IDs...")
  merged$study <- sub("_.*", "", colnames(merged))
  message("  Studies found: ", paste(unique(merged$study), collapse = ", "))
  print(table(merged$study))
} else {
  message("  'study' column exists: ", paste(unique(merged$study), collapse = ", "))
}

merged <- RunHarmony(
	merged,
	group.by.vars = "study",
	verbose       = TRUE
	)

saveRDS(
  merged,
  "results/merged_harmony.rds"
)
message("  Saved: merged_harmony.rds")

# =============================================================
# Post-Harmony: Neighbors + UMAP (use harmony reduction)
# =============================================================
library(Seurat)

options(future.globals.maxSize = 8000 * 1024^2)

message("\nloading Harmony integration...")
merged <- readRDS("results/merged_harmony.rds")

# Filter cells with percent.mt > 10
merged <- subset(
  merged,
  subset = percent_mt <= 10
)



message("Finding neighbors...")
merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims      = 1:30,
  verbose   = FALSE
)

message("Running UMAP...")
merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims      = 1:30,
  verbose   = FALSE
)


# =============================================================
# Clustering at multiple resolutions
# =============================================================
message("Clustering...")
for (res in c(0.1, 0.2, 0.5, 0.8)) {
  merged <- FindClusters(merged, resolution = res, verbose = FALSE)
  message("  Done resolution: ", res)
}

saveRDS(merged, "results/merged_harmony_clusters.rds")
message("  Saved: merged_harmony_clusters.rds")


message("\nDone.")

# =============================================================
# QC plots
# =============================================================
message("Saving UMAP plots...")
dir.create("results/umap_plots", recursive = TRUE, showWarnings = FALSE)

plot_and_save <- function(p, filename, w = 10, h = 8) {
  ggsave(
    file.path("results/umap_plots", filename),
    p, width = w, height = h, dpi = 150
  )
}

# Clustering resolutions
plot_and_save(DimPlot(merged, group.by = "RNA_snn_res.0.1", label = TRUE) + ggtitle("UMAP res 0.1"), "01_int_umap_res0.1.png")
plot_and_save(DimPlot(merged, group.by = "RNA_snn_res.0.2", label = TRUE) + ggtitle("UMAP res 0.2"), "02_int_umap_res0.2.png")
plot_and_save(DimPlot(merged, group.by = "RNA_snn_res.0.5", label = TRUE) + ggtitle("UMAP res 0.5"), "03_int_umap_res0.5.png")
plot_and_save(DimPlot(merged, group.by = "RNA_snn_res.0.8", label = TRUE) + ggtitle("UMAP res 0.8"), "04_int_umap_res0.8.png")

# Study mixing — integration QC
plot_and_save(
  DimPlot(merged, group.by = "study", pt.size = 0.3) + ggtitle("UMAP by Study"),
  "05_umap_study.png", w = 10
)

message("  Saved UMAP plots.")
message("\nDone.")
