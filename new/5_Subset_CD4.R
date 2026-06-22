library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

plot_dir <- "./plots_final_annotations"
dir.create(plot_dir, showWarnings = FALSE)

# ── load final object ─────────────────────────────────────────────
message("loading final object...")
seurat_clean <- readRDS("./main_object_annotated_reclustered.rds")

# ── subset CD4+ effector memory T cells ──────────────────────────
message("subsetting CD4+ effector memory T cells...")
cd4 <- subset(seurat_clean, subset = cell_type == "CD4+ effector memory T cell")
message("CD4+ cells: ", ncol(cd4))

rm(seurat_clean); gc()

# ── recluster using existing harmony reduction ────────────────────
message("available reductions: ", paste(names(cd4@reductions), collapse = ", "))

harmony_reduction <- "harmony"  # adjust if needed
n_dims <- ncol(cd4@reductions[[harmony_reduction]])
message("dims available: ", n_dims)
dims_to_use <- 1:min(20, n_dims)

message("finding neighbors...")
cd4 <- FindNeighbors(cd4, reduction = harmony_reduction, dims = dims_to_use)

# try multiple resolutions
for (res in c(0.1, 0.2, 0.3, 0.5)) {
  cd4 <- FindClusters(cd4, resolution = res,
                       graph.name = "RNA_snn")
  message("res ", res, ": ", length(unique(cd4@meta.data[[paste0("RNA_snn_res.", res)]])), " clusters")
  print(table(cd4@meta.data[[paste0("RNA_snn_res.", res)]]))
}

message("running UMAP...")
cd4 <- RunUMAP(cd4, reduction = harmony_reduction, dims = dims_to_use)

# ── check Tfh / CD4 subtype markers ──────────────────────────────
tfh_markers <- c(
  # pan Tfh
  "CXCR5", "ICOS", "BCL6", "PDCD1", "CD200",
  "TIGIT", "TOX", "TOX2", "MAF",
  # Tfh1
  "CXCR3", "TBX21", "IFNG", "CCR1",
  # Tfh2
  "CCR4", "GATA3", "IL4", "IL13", "ST2",
  # Tfh17
  "CCR6", "RORC", "IL17A", "IL17F",
  # Th1 non-Tfh
  "CCR5", "GZMK", "HAVCR2",
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "IKZF2",
  # naive / memory
  "CCR7", "SELL", "TCF7", "LEF1", "IL7R",
  # activation
  "CD69", "CD44", "HLA-DRA",
  # cytotoxic — should be absent in CD4
  "CD8A", "GZMB", "PRF1"
)

valid_markers <- tfh_markers[tfh_markers %in% rownames(cd4)]
missing       <- tfh_markers[!tfh_markers %in% rownames(cd4)]
if (length(missing) > 0) message("markers not found: ", paste(missing, collapse = ", "))

# ── plots at each resolution ──────────────────────────────────────
for (res in c(0.1, 0.2, 0.3, 0.5)) {

  cluster_col <- paste0("RNA_snn_res.", res)
  n           <- length(unique(cd4@meta.data[[cluster_col]]))

  p_umap <- DimPlot(
    cd4,
    reduction  = "umap",
    group.by   = cluster_col,
    label      = TRUE,
    label.size = 4,
    label.box  = TRUE,
    repel      = TRUE,
    pt.size    = 0.5,
    alpha      = 0.7,
    raster     = FALSE
  ) +
    ggtitle(paste0("CD4+ T cells — res ", res, " (", n, " clusters)")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_dot <- DotPlot(
    cd4,
    features = valid_markers,
    group.by = cluster_col,
    scale    = TRUE
  ) +
    coord_flip() +
    ggtitle(paste0("Tfh markers — res ", res)) +
    theme_classic(base_size = 10) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  combined <- p_umap + p_dot +
    plot_annotation(
      title = paste0("CD4+ T cell subtyping — resolution ", res),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )

  ggsave(
    file.path(plot_dir, paste0("cd4_subtype_res", res, ".pdf")),
    combined, width = 18, height = 8
  )
  message("saved cd4_subtype_res", res, ".pdf")
}

# ── compute markers at res 0.3 to start ──────────────────────────
message("\ncomputing markers at res 0.3...")
Idents(cd4) <- "RNA_snn_res.0.3"

cd4_markers <- FindAllMarkers(
  cd4,
  only.pos   = TRUE,
  min.pct    = 0.25,
  logfc.threshold = 0.5,
  test.use   = "wilcox"
)

cd4_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = pct.1, n = 10) %>%
  print(n = Inf)

write.csv(cd4_markers,
          file.path(plot_dir, "markers_CD4_res0.3.csv"),
          row.names = FALSE)
message("saved markers_CD4_res0.3.csv")

# ── feature plots for key Tfh markers ────────────────────────────
key_features <- c("CXCR5", "BCL6", "ICOS", "CXCR3",
                   "CCR4", "CCR6", "FOXP3", "PDCD1")
key_features <- key_features[key_features %in% rownames(cd4)]

p_feat <- FeaturePlot(
  cd4,
  features   = key_features,
  reduction  = "umap",
  ncol       = 4,
  pt.size    = 0.3,
  order      = TRUE,
  cols       = c("lightgrey", "#E63946")
) &
  theme_classic(base_size = 10) &
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(plot_dir, "cd4_featureplot_tfh_markers.pdf"),
  p_feat, width = 20, height = 10
)
message("saved cd4_featureplot_tfh_markers.pdf")

# ── save CD4 subsetted object ─────────────────────────────────────
saveRDS(cd4, "./cd4_subsetted.rds")
message("saved cd4_subsetted.rds")

message("\n── all done ──")
