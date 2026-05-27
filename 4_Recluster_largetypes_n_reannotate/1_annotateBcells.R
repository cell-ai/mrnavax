setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/results/B_cells/")
library(Seurat)

# Load B-cell object
B_cells <- readRDS("B_cells_seurat.rds")

# Make sure identities are set to resolution 0.3 clusters
Idents(B_cells) <- "RNA_snn_res.0.3"   # change if your column name differs

# Check cluster sizes
table(Idents(B_cells))

# Define clusters to remove
bad.clusters <- c("4", "5", "7", "8", "10")

# Remove bad clusters
B_cells.clean <- subset(B_cells, idents = bad.clusters, invert = TRUE)

# Drop unused factor levels
B_cells.clean$RNA_snn_res.0.3 <- droplevels(B_cells.clean$RNA_snn_res.0.3)

# Re-run clustering workflow
B_cells.clean <- NormalizeData(B_cells.clean)
B_cells.clean <- FindVariableFeatures(B_cells.clean)
B_cells.clean <- ScaleData(B_cells.clean)
B_cells.clean <- RunPCA(B_cells.clean)

B_cells.clean <- FindNeighbors(B_cells.clean, dims = 1:30)
B_cells.clean <- FindClusters(B_cells.clean, resolution = 0.1)
B_cells.clean <- RunUMAP(B_cells.clean, dims = 1:30)

DimPlot(B_cells.clean, reduction = "umap")

pan_b <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1")
pan_t <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2")

# Plot individual markers
FeaturePlot(
  B_cells.clean,
  features = c(pan_b, pan_t),
  reduction = "umap",
  ncol = 4
)

# DotPlot by cluster/annotation
DotPlot(
  B_cells.clean,
  features = list(
    Pan_B = pan_b,
    Pan_T = pan_t
  )
) + RotatedAxis()

# Find markers after reclustering
B_cells.clean.markers <- FindAllMarkers(
  B_cells.clean,
  only.pos = TRUE,
  min.pct = 0.50,
  logfc.threshold = 0.25
)

top10 <- B_cells.clean.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster)

write.csv(top10, "B_cells_clean_reclustered_markers.csv")


# verify
table(Idents(B_cells.clean))

############################Reclustering again
# Define clusters to remove
bad.clusters <- c("6")

# Remove bad clusters
B_cells.clean <- subset(B_cells.clean, idents = bad.clusters, invert = TRUE)

# Drop unused factor levels
B_cells.clean$RNA_snn_res.0.1 <- droplevels(B_cells.clean$RNA_snn_res.0.1)

# Re-run clustering workflow
B_cells.clean <- NormalizeData(B_cells.clean)
B_cells.clean <- FindVariableFeatures(B_cells.clean)
B_cells.clean <- ScaleData(B_cells.clean)
B_cells.clean <- RunPCA(B_cells.clean)

B_cells.clean <- FindNeighbors(B_cells.clean, dims = 1:30)
B_cells.clean <- FindClusters(B_cells.clean, resolution = 0.1)
B_cells.clean <- RunUMAP(B_cells.clean, dims = 1:30)

DimPlot(B_cells.clean, reduction = "umap")

pan_b <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1")
pan_t <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2")

# Plot individual markers
FeaturePlot(
  B_cells.clean,
  features = c(pan_b, pan_t),
  reduction = "umap",
  ncol = 4
)

# DotPlot by cluster/annotation
DotPlot(
  B_cells.clean,
  features = list(
    Pan_B = pan_b,
    Pan_T = pan_t
  )
) + RotatedAxis()

# Find markers after reclustering
B_cells.clean.markers <- FindAllMarkers(
  B_cells.clean,
  only.pos = TRUE,
  min.pct = 0.50,
  logfc.threshold = 0.25
)

top10 <- B_cells.clean.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster)

write.csv(top10, "B_cells_clean_reclustered_markers.csv")

# ── 3. Plot ───────────────────────────────────────────────────────────────
DimPlot(
  B_cells.clean,
  reduction  = "umap",
  label      = TRUE,
  label.size = 3.5,
  repel      = TRUE,
  pt.size    = 0.4
) +
  ggtitle("B cell subclusters") +
  theme(legend.text = element_text(size = 9))


# Save outputs
saveRDS(B_cells.clean, "B_cells_clean_reclustered_res0.1.rds")
write.csv(
  B_cells.clean.markers,
  "B_cells_clean_reclustered_res0.3_markers.csv",
  row.names = FALSE
)

pan_b <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1")
pan_t <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2")

# Plot individual markers
FeaturePlot(
  B_cells.clean,
  features = c(pan_b, pan_t),
  reduction = "umap",
  ncol = 4
)

# DotPlot by cluster/annotation
DotPlot(
  B_cells.clean,
  features = list(
    Pan_B = pan_b,
    Pan_T = pan_t
  )
) + RotatedAxis()
