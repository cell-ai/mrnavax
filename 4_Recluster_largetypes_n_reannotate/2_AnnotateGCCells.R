setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/results/GCB_Plasma/")
library(Seurat)

# Load T-cell object
GC_cells <- readRDS("GCB_Plasma_seurat.rds")

# Make sure identities are set to resolution 0.3 clusters
Idents(GC_cells) <- "RNA_snn_res.0.3"   # change if your column name differs

# Check cluster sizes
table(Idents(GC_cells))

# Define clusters to remove
bad.clusters <- c("5", "7", "9") #Stromal, t-cell, possible contaminant

# Remove bad clusters
GC_cells.clean <- subset(GC_cells, idents = bad.clusters, invert = TRUE)

# Drop unused factor levels
GC_cells.clean$RNA_snn_res.0.3 <- droplevels(GC_cells.clean$RNA_snn_res.0.3)

# Re-run clustering workflow
GC_cells.clean <- NormalizeData(GC_cells.clean)
GC_cells.clean <- FindVariableFeatures(GC_cells.clean)
GC_cells.clean <- ScaleData(GC_cells.clean)
GC_cells.clean <- RunPCA(GC_cells.clean)

GC_cells.clean <- FindNeighbors(GC_cells.clean, dims = 1:30)
GC_cells.clean <- FindClusters(GC_cells.clean, resolution = 0.1)
GC_cells.clean <- RunUMAP(GC_cells.clean, dims = 1:30)

# Find markers after reclustering
GC_cells.clean.markers <- FindAllMarkers(
  GC_cells.clean,
  only.pos = TRUE,
  min.pct = 0.50,
  logfc.threshold = 0.25
)

top10 <- GC_cells.clean.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster)

write.csv(top10, "GC_cells_clean_reclustered_markers.csv")

table(Idents(GC_cells.clean))

pan_b <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1")
pan_t <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2")

# Plot individual markers
FeaturePlot(
  GC_cells.clean,
  features = c(pan_b, pan_t),
  reduction = "umap",
  ncol = 4
)

# DotPlot by cluster/annotation
DotPlot(
  GC_cells.clean,
  features = list(
    Pan_B = pan_b,
    Pan_T = pan_t
  )
) + RotatedAxis()


# ── 1. Define cluster labels ──────────────────────────────────────────────
cluster_labels <- c(
  "0"  = "Naive/Stem-like T",
  "1"  = "CD161+ Memory T",
  "2"  = "Cytotoxic CD8 T",
  "3"  = "Activated Cytotoxic T/NK-like",
  "4"  = "Treg",
  "5"  = "Tfh/Tph-like Vaccine-responsive T",
  "6"  = "Activated/Effector Treg",
  "7"  = "IFN-stimulated T",
  "8"  = "Stress-response/Transitional T",
  "9"  = "IL10+ Regulatory T",
  "10" = "Gamma-delta T",
  "11" = "Ig-high T",
  "12" = "MHC-II+/CD74+ Activated T"
)

# ── 2. Map labels to Seurat object ────────────────────────────────────────
cell_identity <- cluster_labels[
  as.character(T_cells.clean$seurat_clusters)
]

names(cell_identity) <- colnames(T_cells.clean)

T_cells.clean <- AddMetaData(
  T_cells.clean,
  metadata = cell_identity,
  col.name = "cell_identity"
)

Idents(T_cells.clean) <- "cell_identity"

# verify
table(Idents(T_cells.clean))

# ── 3. Plot ───────────────────────────────────────────────────────────────
DimPlot(
  T_cells.clean,
  reduction  = "umap",
  label      = TRUE,
  label.size = 3.5,
  repel      = TRUE,
  pt.size    = 0.4
) +
  ggtitle("T cell subclusters") +
  theme(
    legend.text = element_text(size = 9)
  )

# ── 4. Save outputs ───────────────────────────────────────────────────────
saveRDS(
  T_cells.clean,
  "T_cells_clean_reclustered_res0.3_annotated.rds"
)

write.csv(
  T_cells.clean.markers,
  "T_cells_clean_reclustered_res0.3_markers.csv",
  row.names = FALSE
)

pan_b <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1", "CD79A")
pan_t <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2")

# Plot individual markers
FeaturePlot(
  T_cells.clean,
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
