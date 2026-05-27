setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/results/Proliferating/")
library(Seurat)

# Load T-cell object
Pr_cells <- readRDS("Proliferating_seurat.rds")

# Make sure identities are set to resolution 0.3 clusters
Idents(Pr_cells) <- "RNA_snn_res.0.2"   # change if your column name differs

# Check cluster sizes
table(Idents(Pr_cells))

# Define clusters to remove
bad.clusters <- c("3", "4", "6", "5", "8") #contaminant

# Remove bad clusters
Pr_cells.clean <- subset(Pr_cells, idents = bad.clusters, invert = TRUE)

# Drop unused factor levels
Pr_cells.clean$RNA_snn_res.0.2 <- droplevels(Pr_cells.clean$RNA_snn_res.0.2)

# Re-run clustering workflow
Pr_cells.clean <- NormalizeData(Pr_cells.clean)
Pr_cells.clean <- FindVariableFeatures(Pr_cells.clean)
Pr_cells.clean <- ScaleData(Pr_cells.clean)
Pr_cells.clean <- RunPCA(Pr_cells.clean)

Pr_cells.clean <- FindNeighbors(Pr_cells.clean, dims = 1:30)
Pr_cells.clean <- FindClusters(Pr_cells.clean, resolution = 0.1)
Pr_cells.clean <- RunUMAP(Pr_cells.clean, dims = 1:30)

# Find markers after reclustering
Pr_cells.clean.markers <- FindAllMarkers(
  Pr_cells.clean,
  only.pos = TRUE,
  min.pct = 0.50,
  logfc.threshold = 0.25
)

top10 <- Pr_cells.clean.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster)

write.csv(top10, "Pr_cells_clean_reclustered_markers.csv")

table(Idents(Pr_cells.clean))

pan_b <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1")
pan_t <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2")

# Plot individual markers
FeaturePlot(
  Pr_cells.clean,
  features = c(pan_b, pan_t),
  reduction = "umap",
  ncol = 4
)

# DotPlot by cluster/annotation
DotPlot(
  Pr_cells.clean,
  features = list(
    Pan_B = pan_b,
    Pan_T = pan_t
  )
) + RotatedAxis()

