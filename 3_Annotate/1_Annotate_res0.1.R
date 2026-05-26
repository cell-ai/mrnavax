library(Seurat)
library(dplyr)
library(ggplot2)

# =============================================================
# Load RDS
# =============================================================
message("Loading merged_harmony_clusters.rds...")
merged <- readRDS("merged_harmony_clusters.rds")
message("  Loaded.")

# =============================================================
# Join layers (Seurat v5 requirement)
# =============================================================
message("Joining layers...")
merged <- JoinLayers(merged)

# =============================================================
# Find all markers at res 0.1
# =============================================================
Idents(merged) <- "RNA_snn_res.0.1"
message("Running FindAllMarkers at res 0.1...")

markers_all <- FindAllMarkers(
  merged,
  only.pos            = TRUE,
  logfc.threshold     = 0.5,
  min.pct             = 0.25,
  max.cells.per.ident = 500
)

write.csv(markers_all, "results/markers_res0.1.csv", row.names = FALSE)
message("  Saved: markers_res0.1.csv")

# =============================================================
# Top 5 markers per cluster
# =============================================================
top5 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

write.csv(top5, "results/top5_markers_res0.1.csv", row.names = FALSE)
message("  Saved: top5_markers_res0.1.csv")

# =============================================================
# Top 10 markers per cluster
# =============================================================
top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10, "results/top10_markers_res0.1.csv", row.names = FALSE)
message("  Saved: top10_markers_res0.1.csv")


# =============================================================
# Bubble plot
# =============================================================
message("Saving bubble plot...")
dir.create("results/umap_plots", recursive = TRUE, showWarnings = FALSE)

top_genes <- unique(top5$gene)

png(
  "results/umap_plots/bubble_top5_markers.png",
  width  = 1800,
  height = 800
)
DotPlot(
  merged,
  features = top_genes,
  dot.scale = 6
) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Top 5 markers per cluster (res 0.1)")
dev.off()

message("  Saved: bubble_top5_markers.png")
message("\nDone.")
