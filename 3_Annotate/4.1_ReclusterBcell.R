#!/usr/bin/env Rscript

setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

# ── 1. Load RDS ──────────────────────────────────────────────
message("Carregando objeto...")
seurat_obj <- readRDS("merged_harmony_clusters.rds")

# ── 2. Set active identity ────────────────────────────────────
cluster_col <- "RNA_snn_res.0.1"
Idents(seurat_obj) <- cluster_col

# ── 3. Cluster annotation ─────────────────────────────────────
cluster_annotation <- list(
  "B_cells" = c(1, 2, 7, 11)   # B maduro, B naive, B transicional, B memoria
)

# ── 4. Funcao de reclusterizacao ──────────────────────────────
recluster_subset <- function(seurat_obj,
                             cluster_col,
                             clusters_to_keep,
                             subset_name,
                             resolution      = 0.3,
                             n_pcs           = 15,
                             n_variable      = 2000,
                             min_pct         = 0.1,
                             logfc_threshold = 0.25) {

  message("\n========== ", subset_name, " ==========")

  # 4.1 Subsetar
  cells_keep <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data[[cluster_col]] %in% clusters_to_keep
  ]
  sub <- subset(seurat_obj, cells = cells_keep)

  # Seurat v5: juntar layers separados por amostra
  sub <- JoinLayers(sub)
  message("Celulas no subset: ", ncol(sub))

  # 4.2 Renormalizar e reescalar
  sub <- NormalizeData(sub, verbose = FALSE)
  message("Data Normalized...")
  sub <- FindVariableFeatures(sub, nfeatures = n_variable, verbose = FALSE)
  message("Found Variable Features...")
  sub <- ScaleData(sub, verbose = FALSE)
  message("Data scaled...")
  # 4.3 PCA + re-Harmony
  sub <- RunPCA(sub, npcs = n_pcs, verbose = FALSE)
  message("PCA done...") 
  sub <- RunHarmony(sub, "orig.ident", reduction.use = "pca", dims.use = 1:n_pcs)
  message("Harmony done...")
  # 4.4 Clustering
  sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  message("Found neighbors...")
  sub <- FindClusters(sub, resolution = resolution, verbose = FALSE)
  message("Found clusters...")
  # 4.5 UMAP
  message("Rodando UMAP em ", ncol(sub), " celulas...")
  sub <- RunUMAP(sub,
                 reduction   = "harmony",
                 dims        = 1:n_pcs,
                 umap.method = "uwot",
                 n.neighbors = 30L,
                 n.epochs    = 200L,
                 approx.pow  = TRUE,
                 verbose     = FALSE)

  # 4.6 FindAllMarkers
  message("Rodando FindAllMarkers...")
  markers <- FindAllMarkers(
    sub,
    only.pos        = TRUE,
    min.pct         = min_pct,
    logfc.threshold = logfc_threshold,
    test.use        = "wilcox"
  )

  markers <- markers %>%
    filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))

  # 4.7 Salvar resultados
  out_dir <- paste0("results/", subset_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(sub, file = paste0(out_dir, "/", subset_name, "_seurat.rds"))
  write.csv(markers, file = paste0(out_dir, "/", subset_name, "_markers.csv"),
            row.names = FALSE)

  # 4.8 Plot UMAP
  png(paste0(out_dir, "/", subset_name, "_umap.png"), width = 1600, height = 700, res = 150)
  p <- DimPlot(sub, reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggtitle(paste0(subset_name, " - res ", resolution))
  print(p)
  dev.off()

  message("Salvos em: ", out_dir)
  return(list(seurat = sub, markers = markers))
}

# ── 5. Rodar B cells ──────────────────────────────────────────
b_cells <- recluster_subset(
  seurat_obj       = seurat_obj,
  cluster_col      = cluster_col,
  clusters_to_keep = cluster_annotation$B_cells,
  subset_name      = "B_cells",
  resolution       = 0.3,
  n_pcs            = 15
)

message("\nConcluido! Resultados em results/B_cells/")
