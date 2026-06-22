setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(Seurat)

plot_dir <- "plots_final_annotations/"

adata <- readRDS("cd4_subsetted.rds")
colnames(adata[[]])

DimPlot(adata, reduction = "umap", group.by = "RNA_snn_res.0.8")
# ── quick look at cluster sizes ───────────────────────────────────
message("cluster sizes at res 0.1:")
print(table(adata$RNA_snn_res.0.1))

message("\ncluster sizes at res 0.3:")
print(table(adata$RNA_snn_res.0.3))

# ══════════════════════════════════════════════════════════════════
# MARKERS AT RESOLUTION 0.1
# ══════════════════════════════════════════════════════════════════
message("\ncomputing markers at res 0.1...")
Idents(adata) <- "RNA_snn_res.0.1"

markers_01 <- FindAllMarkers(
  adata,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5,
  test.use        = "wilcox"
)

# top 10 by pct.1
message("\nTop 10 markers per cluster at res 0.1 (by pct.1):")
markers_01 %>%
  group_by(cluster) %>%
  slice_max(order_by = pct.1, n = 10) %>%
  print(n = Inf)

write.csv(markers_01,
          file.path(plot_dir, "markers_CD4_res0.1.csv"),
          row.names = FALSE)
message("saved markers_CD4_res0.1.csv")

# ══════════════════════════════════════════════════════════════════
# MARKERS AT RESOLUTION 0.3
# ══════════════════════════════════════════════════════════════════
message("\ncomputing markers at res 0.3...")
Idents(adata) <- "RNA_snn_res.0.3"

markers_03 <- FindAllMarkers(
  adata,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5,
  test.use        = "wilcox"
)

message("\nTop 10 markers per cluster at res 0.3 (by pct.1):")
markers_03 %>%
  group_by(cluster) %>%
  slice_max(order_by = pct.1, n = 10) %>%
  print(n = Inf)

plot_dir <- "plots_final_annotations/"

write.csv(markers_03,
          file.path(plot_dir, "markers_CD4_res0.3.csv"),
          row.names = FALSE)
message("saved markers_CD4_res0.3.csv")

############################################################################

cd4_c2 <- subset(adata, subset = RNA_snn_res.0.3 == "2")

#----------------------------------
# 2. Re-cluster
#----------------------------------
DefaultAssay(cd4_c2) <- "RNA"

cd4_c2 <- FindVariableFeatures(cd4_c2)
cd4_c2 <- ScaleData(cd4_c2)
cd4_c2 <- RunPCA(cd4_c2)

cd4_c2 <- FindNeighbors(cd4_c2, dims = 1:30)

cd4_c2 <- FindClusters(cd4_c2, resolution = c(0.1, 0.3))

cd4_c2 <- RunUMAP(cd4_c2, dims = 1:30)

# check cluster numbers
DimPlot(cd4_c2, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(cd4_c2, group.by = "RNA_snn_res.0.3", label = TRUE)

Th1_sig <-list(c("TBX21", "ID2", "IFNG", "IL2RB", "STAT1", "CXCR3", "CCR5"))

#edit here stat5 to Stat5a and Stat5b
Th2_sig <- list(c("GATA3", "IRF4", "STAT5a", "STAT5b", "IL1RL1", "CXCR4", "IL4"))

#edit here Il6r to Il6ra
Th17_sig <- list(c("IL17A", "IL17F", "RORC", "RORA", "CCR6", "IL23R", "IL6RA", "STAT3"))


signatures<- c(Th1_sig, Th2_sig, Th17_sig)

cd4_c2 <- AddModuleScore_UCell(
  cd4_c2,
  features = signatures
)

scores_01 <- cd4_c2@meta.data %>%
  mutate(cluster = cd4_c2$RNA_snn_res.0.1) %>%
  group_by(cluster) %>%
  summarise(
    Sig1 = mean(signature_1_UCell),
    Sig2 = mean(signature_2_UCell),
    Sig3 = mean(signature_3_UCell)
  )

mat01 <- as.data.frame(scores_01)
rownames(mat01) <- mat01$cluster
mat01$cluster <- NULL

library(pheatmap)

pheatmap(
  as.matrix(mat01),
  scale = "column",
  main = "Resolution 0.1"
)

scores_03 <- cd4_c2@meta.data %>%
  mutate(cluster = cd4_c2$RNA_snn_res.0.3) %>%
  group_by(cluster) %>%
  summarise(
    Sig1 = mean(signature_1_UCell),
    Sig2 = mean(signature_2_UCell),
    Sig3 = mean(signature_3_UCell)
  )

mat03 <- as.data.frame(scores_03)
rownames(mat03) <- mat03$cluster
mat03$cluster <- NULL

pheatmap(
  as.matrix(mat03),
  scale = "column",
  main = "Resolution 0.3"
)


p <- FeaturePlot(
  cd4_c2,
  features = "signature_1_UCell",
  reduction = "umap",
  order = TRUE
)

LabelClusters(p, id = "ident")


library(UCell)

tfh_signatures <- list(
  TFH1  = c("CXCR3", "TBX21", "IFNG", "CCL5", "CTSW", "GZMM"),
  TFH2  = c("GATA3", "CCR4", "CXCR4", "MAF", "IL4", "IL13"),
  TFH17 = c("RORC", "CCR6", "RORA", "IL17A", "IL21", "KLRB1"),
  Tfr   = c("FOXP3", "IKZF2", "CTLA4", "IL2RA", "TIGIT", "FCRL3")
)

seurat_obj <- AddModuleScore_UCell(seurat_obj, features = tfh_signatures)

score_cols <- c("TFH1_UCell", "TFH2_UCell", "TFH17_UCell", "Tfr_UCell")

ucell_summary <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells        = n(),
    TFH1_mean      = round(mean(TFH1_UCell),  4),
    TFH2_mean      = round(mean(TFH2_UCell),  4),
    TFH17_mean     = round(mean(TFH17_UCell), 4),
    Tfr_mean       = round(mean(Tfr_UCell),   4),
    .groups = "drop"
  ) %>%
  mutate(top_identity = apply(.[, c("TFH1_mean","TFH2_mean","TFH17_mean","Tfr_mean")], 1,
                              function(x) c("TFH1","TFH2","TFH17","Tfr")[which.max(x)]))

print(ucell_summary)

cluster_labels <- c(
  "0" = "Resting/Memory TFH",
  "1" = "Activated TFH (OX40+, S100+)",
  "2" = "TFH1 (CXCR3+, CCL5+, NELL2+)",
  "3" = "Tfr (follicular Treg)",
  "4" = "IFN-stimulated TFH (ISG-hi)"
)

# If still failing, do it via index
seurat_obj$tfh_subtype <- plyr::mapvalues(
  x    = as.character(seurat_obj$RNA_snn_res.0.3),
  from = names(cluster_labels),
  to   = cluster_labels
)
