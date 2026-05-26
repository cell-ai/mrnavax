#!/usr/bin/env Rscript

setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

message("Carregando objeto...")
seurat_obj <- readRDS("merged_harmony_clusters.rds")

cluster_col <- "RNA_snn_res.0.1"
Idents(seurat_obj) <- cluster_col

cluster_annotation <- list(
  "T_cells" = c(0, 6, 10)   # T naive, Tfh, T memoria
)

recluster_subset <- function(seurat_obj,
                             cluster_col,
                             clusters_to_keep,
                             subset_name,
                             resolutions     = c(0.1, 0.3, 0.5),
                             n_pcs           = 20,
                             n_variable      = 2000,
                             min_pct         = 0.1,
                             logfc_threshold = 0.25) {

  message("\n========== ", subset_name, " ==========")

  out_dir <- paste0("results/", subset_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cells_keep <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data[[cluster_col]] %in% clusters_to_keep
  ]

  sub <- subset(seurat_obj, cells = cells_keep)

  sub <- JoinLayers(sub)
  DefaultAssay(sub) <- "RNA"

  message("Celulas no subset: ", ncol(sub))

  sub <- NormalizeData(sub, verbose = FALSE)
  message("Data normalized...")

  sub <- FindVariableFeatures(sub, nfeatures = n_variable, verbose = FALSE)
  message("Found variable features...")

  sub <- ScaleData(sub, verbose = FALSE)
  message("Data scaled...")

  sub <- RunPCA(sub, npcs = n_pcs, verbose = FALSE)
  message("PCA done...")

  sub <- RunHarmony(
    sub,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    dims.use      = 1:n_pcs,
    verbose       = FALSE
  )
  message("Harmony done...")

  sub <- FindNeighbors(
    sub,
    reduction = "harmony",
    dims      = 1:n_pcs,
    verbose   = FALSE
  )
  message("Found neighbors...")

  results <- list()

  for (res in resolutions) {

    res_label <- gsub("\\.", "_", as.character(res))
    message("\n--- Rodando resolucao ", res, " ---")

    sub_res <- FindClusters(
      sub,
      resolution = res,
      verbose    = FALSE
    )

    Idents(sub_res) <- "seurat_clusters"
    message("Found clusters...")

    message("Rodando UMAP em ", ncol(sub_res), " celulas...")
    sub_res <- RunUMAP(
      sub_res,
      reduction   = "harmony",
      dims        = 1:n_pcs,
      umap.method = "uwot",
      n.neighbors = 30L,
      n.epochs    = 200L,
      approx.pow  = TRUE,
      verbose     = FALSE
    )

    message("Rodando FindAllMarkers...")
    markers <- FindAllMarkers(
      sub_res,
      only.pos        = TRUE,
      min.pct         = min_pct,
      logfc.threshold = logfc_threshold,
      test.use        = "wilcox"
    )

    markers <- markers %>%
      filter(p_val_adj < 0.05) %>%
      arrange(cluster, desc(avg_log2FC))

    saveRDS(
      sub_res,
      file = paste0(out_dir, "/", subset_name, "_res", res_label, "_seurat.rds")
    )

    write.csv(
      markers,
      file = paste0(out_dir, "/", subset_name, "_res", res_label, "_markers.csv"),
      row.names = FALSE
    )

    png(
      paste0(out_dir, "/", subset_name, "_res", res_label, "_umap.png"),
      width = 1600,
      height = 700,
      res = 150
    )

    p <- DimPlot(
      sub_res,
      reduction = "umap",
      label = TRUE,
      pt.size = 0.5
    ) +
      ggtitle(paste0(subset_name, " - res ", res))

    print(p)
    dev.off()

    results[[paste0("res_", res_label)]] <- list(
      seurat = sub_res,
      markers = markers
    )

    message("Salvos resultados da resolucao ", res, " em: ", out_dir)
  }

  return(results)
}

t_cells <- recluster_subset(
  seurat_obj       = seurat_obj,
  cluster_col      = cluster_col,
  clusters_to_keep = cluster_annotation$T_cells,
  subset_name      = "T_cells",
  resolutions      = c(0.1, 0.3, 0.5),
  n_pcs            = 20
)

message("\nConcluido! Resultados em results/T_cells/")
