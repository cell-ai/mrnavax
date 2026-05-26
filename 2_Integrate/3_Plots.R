library(Seurat)
library(ggplot2)
library(pheatmap)
# =============================================================
# Load RDS
# =============================================================
message("Loading merged_harmony_clusters.rds...")
merged <- readRDS("results/merged_harmony_clusters.rds")
message("  Loaded.")

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
plot_and_save(DimPlot(merged, group.by = "sample", label = TRUE) + ggtitle("Sample"), "05_int_umap_sample.png")

# Study mixing — integration QC
plot_and_save(
  DimPlot(merged, group.by = "study", pt.size = 0.3) + ggtitle("UMAP by Study"),
  "05_umap_study.png", w = 10
)

message("  Saved UMAP plots.")
message("\nDone.")

# Calculate ribosomal percentage if needed
merged[["percent.ribo"]] <- PercentageFeatureSet(
    merged,
    pattern = "^RPL|^RPS"
)

library(ggplot2)

# ----------------------------
# Plot 1: MT vs Ribosomal
# ----------------------------
p1 <- FeatureScatter(
    merged,
    feature1 = "percent_mt",
    feature2 = "percent.ribo",
    group.by = "seurat_clusters"
)

ggsave(
    "MT_vs_Ribosomal_scatter.png",
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
)

# ----------------------------
# Plot 2: MT vs nFeature
# ----------------------------
p2 <- FeatureScatter(
    merged,
    feature1 = "percent_mt",
    feature2 = "nFeature_RNA",
    group.by = "seurat_clusters"
)

ggsave(
    "MT_vs_nFeature_scatter.png",
    plot = p2,
    width = 8,
    height = 6,
    dpi = 300
)

# ----------------------------
# Plot 3: Ribosomal vs nFeature
# ----------------------------
p3 <- FeatureScatter(
    merged,
    feature1 = "percent.ribo",
    feature2 = "nFeature_RNA",
    group.by = "seurat_clusters"
)

ggsave(
    "Ribosomal_vs_nFeature_scatter.png",
    plot = p3,
    width = 8,
    height = 6,
    dpi = 300
)

# ----------------------------
# Optional violin plot for QC
# ----------------------------
p4 <- VlnPlot(
    merged,
    features=c(
        "percent_mt",
        "percent.ribo",
        "nFeature_RNA",
        "nCount_RNA"
    ),
    group.by="seurat_clusters",
    pt.size=0
)

ggsave(
    "QC_violin_clusters.png",
    plot = p4,
    width = 12,
    height = 6,
    dpi = 300
)

# ----------------------------
# 1. Check metadata columns
# ----------------------------
colnames(merged@meta.data)

# Choose sample column
# Change this if your sample column has another name
sample_col <- "orig.ident"

# Check available samples
table(merged@meta.data[[sample_col]])

# ----------------------------
# 2. Count cells per sample per cluster
# ----------------------------
counts_mat <- table(
  merged@meta.data[[sample_col]],
  Idents(merged)
)

counts_mat

write.csv(
  as.data.frame.matrix(counts_mat),
  file = "sample_by_cluster_cell_counts.csv"
)

# ----------------------------
# 3. Proportion of each cluster from each sample
# ----------------------------
prop_mat <- prop.table(counts_mat, margin = 2)

prop_mat

write.csv(
  as.data.frame.matrix(prop_mat),
  file = "sample_by_cluster_proportions.csv"
)

# ----------------------------
# 4. Specifically inspect cluster 0
# ----------------------------
cluster0_counts <- table(
  merged@meta.data[[sample_col]][Idents(merged) == "0"]
)

cluster0_props <- prop.table(cluster0_counts)

cluster0_counts
cluster0_props

write.csv(
  as.data.frame(cluster0_counts),
  file = "cluster0_sample_counts.csv"
)

write.csv(
  as.data.frame(cluster0_props),
  file = "cluster0_sample_proportions.csv"
)

# ----------------------------
# 5. Specifically inspect cluster 14
# ----------------------------
cluster14_counts <- table(
  merged@meta.data[[sample_col]][Idents(merged) == "14"]
)

cluster14_props <- prop.table(cluster14_counts)

cluster14_counts
cluster14_props

write.csv(
  as.data.frame(cluster14_counts),
  file = "cluster14_sample_counts.csv"
)

write.csv(
  as.data.frame(cluster14_props),
  file = "cluster14_sample_proportions.csv"
)

# ----------------------------
# 6. Heatmap: sample contribution to each cluster
# ----------------------------
png(
  filename = "sample_by_cluster_proportion_heatmap.png",
  width = 1800,
  height = 1200,
  res = 200
)

pheatmap(
  prop_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Sample contribution per cluster"
)

dev.off()

# ----------------------------
# 7. Barplot: sample composition by cluster
# ----------------------------
df <- as.data.frame(counts_mat)
colnames(df) <- c("sample", "cluster", "cells")

p <- ggplot(df, aes(x = cluster, y = cells, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  ylab("Proportion of cells") +
  xlab("Cluster") +
  ggtitle("Sample composition by cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "sample_composition_by_cluster_barplot.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)
