#!/usr/bin/env Rscript
# =============================================================
# 1.Processing.R
# Owner: Izabela Mamede
# Date:  2026-05-17
# Altered by Marcela
# =============================================================

setwd("/media/csbl/sandbox-SSD-3/mra_vaccine")

STUDY <- "GSE195673"

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(harmony)
})

# ── Parallel settings ─────────────────────────────────────────

N_THREADS <- 16   # change to your CPU count

future::plan("multicore", workers = N_THREADS)

# increase maximum object size passed between workers
options(future.globals.maxSize = 50 * 1024^3)  # 50 GB

message("Running with ", N_THREADS, " threads")

dir.create("results/qc_plots",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/umap_plots", recursive = TRUE, showWarnings = FALSE)

extract_or_unknown <- function(x, pattern) {
  m <- regmatches(x, regexpr(pattern, x))
  if (length(m) == 0 || is.na(m) || m == "") "unknown" else m
}

split_aggr_into_samples <- function(matrix_path, features_path,
                                    barcodes_path, aggregation_csv) {

  message("Reading matrix...")
  mat      <- as(readMM(matrix_path), "CsparseMatrix")
  features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  agg      <- read.csv(aggregation_csv, stringsAsFactors = FALSE)

  if (!all(c("sample_id", "library_id") %in% colnames(agg))) {
    stop("aggregation_csv must contain columns: sample_id and library_id")
  }

  rownames(mat) <- ifelse(
    !is.na(features$V2) & features$V2 != "",
    features$V2,
    features$V1
  )

  if (any(duplicated(rownames(mat)))) {
    message("  Deduplicating gene names...")
    rownames(mat) <- make.unique(rownames(mat))
  }

  colnames(mat) <- barcodes$V1

  barcode_sample_idx <- suppressWarnings(as.integer(sub(".*-", "", barcodes$V1)))

  if (any(is.na(barcode_sample_idx))) {
    stop("Some barcodes do not end with a numeric sample suffix, e.g. '-1', '-2'.")
  }

  sample_objects <- list()

  for (idx in seq_len(nrow(agg))) {

    sample_id <- agg$sample_id[idx]
    lib       <- agg$library_id[idx]
    cell_mask <- barcode_sample_idx == idx

    if (sum(cell_mask) == 0) {
      message("  Skipping ", sample_id, " — no cells found for index ", idx)
      next
    }

    sub_mat <- mat[, cell_mask, drop = FALSE]

    obj <- CreateSeuratObject(
      counts       = sub_mat,
      project      = sample_id,
      min.cells    = 3,
      min.features = 200
    )

    obj$sample <- sample_id
    obj$gse <- sub("_aggregation.csv", "", basename(aggregation_csv))

    obj$tissue <- dplyr::case_when(
      grepl("_PB_", lib) ~ "PBMC",
      grepl("_LN_", lib) ~ "lymph_node",
      TRUE ~ "unknown"
    )

    obj$timepoint <- extract_or_unknown(lib, "d\\d+")

    obj$donor <- if (grepl("ELAB-WU\\d+", lib)) {
      sub(".*ELAB-(WU\\d+).*", "\\1", lib)
    } else {
      "unknown"
    }

    obj <- RenameCells(obj, add.cell.id = sample_id)

    message("  Sample: ",     sample_id,
            " | Tissue: ",    unique(obj$tissue),
            " | Timepoint: ", unique(obj$timepoint),
            " | Donor: ",     unique(obj$donor),
            " | Cells: ",     ncol(obj))

    sample_objects[[sample_id]] <- obj
  }

  return(sample_objects)
}

seurat_list <- split_aggr_into_samples(
  matrix_path     = paste0("h5ad_geo/", STUDY, "_matrix.mtx.gz"),
  features_path   = paste0("h5ad_geo/", STUDY, "_features.tsv.gz"),
  barcodes_path   = paste0("h5ad_geo/", STUDY, "_barcodes.tsv.gz"),
  aggregation_csv = paste0("h5ad_geo/", STUDY, "_aggregation.csv")
)

message("Total samples loaded: ", length(seurat_list))

message("Merging samples...")

if (length(seurat_list) == 0) {
  stop("No samples were loaded. Check aggregation CSV and barcode suffixes.")
}

if (length(seurat_list) == 1) {
  seurat_merged <- seurat_list[[1]]
} else {
  seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1], project = STUDY)
}

if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])
}

message("Computing QC metrics...")
seurat_merged[["percent_mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")
seurat_merged$log10GenesPerUMI <- log10(seurat_merged$nFeature_RNA) / log10(seurat_merged$nCount_RNA)

message("Plotting QC violins...")
p1 <- VlnPlot(
  seurat_merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
  group.by = "sample",
  pt.size  = 0,
  raster   = TRUE,
  ncol     = 1
) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(
  paste0("results/qc_plots/", STUDY, "_01_violin_qc.png"),
  p1,
  width = 14,
  height = 12,
  dpi = 150
)

summary_stats <- seurat_merged@meta.data %>%
  group_by(sample, gse, tissue, timepoint, donor) %>%
  summarise(
    n_cells         = n(),
    nFeature_median = median(nFeature_RNA),
    nCount_median   = median(nCount_RNA),
    mt_median       = round(median(percent_mt), 3),
    mt_max          = round(max(percent_mt), 3),
    .groups = "drop"
  ) %>%
  arrange(timepoint)

print(summary_stats)

write.csv(
  summary_stats,
  paste0("results/qc_plots/", STUDY, "_summary_stats.csv"),
  row.names = FALSE
)

message("Filtering cells...")
seurat_qc <- subset(
  seurat_merged,
  subset = nCount_RNA       >= 500   &
           nCount_RNA       <= 60000 &
           nFeature_RNA     >= 500   &
           nFeature_RNA     <= 8000  &
           log10GenesPerUMI >= 0.80  &
           percent_mt       <  20
)

message("Cells after QC filter: ", ncol(seurat_qc), " (from ", ncol(seurat_merged), ")")

DefaultAssay(seurat_qc) <- "RNA"

message("Normalizing...")
seurat_qc <- NormalizeData(seurat_qc, verbose = FALSE)
seurat_qc <- FindVariableFeatures(seurat_qc, nfeatures = 2000, verbose = FALSE)
seurat_qc <- ScaleData(seurat_qc, verbose = FALSE)
seurat_qc <- RunPCA(seurat_qc, npcs = 50, verbose = FALSE)

message("Running Harmony...")
seurat_harmony <- RunHarmony(
  seurat_qc,
  group.by.vars = "sample",
  verbose = FALSE
)

message("Running UMAP / tSNE...")
seurat_harmony <- RunUMAP(
  seurat_harmony,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

seurat_harmony <- RunTSNE(
  seurat_harmony,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

seurat_harmony <- FindNeighbors(
  seurat_harmony,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

message("Clustering...")
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.1, verbose = FALSE)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.2, verbose = FALSE)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5, verbose = FALSE)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.8, verbose = FALSE)

message("Saving UMAP plots...")

plot_and_save <- function(p, filename, w = 10, h = 8) {
  ggsave(
    file.path("results/umap_plots", paste0(STUDY, "_", filename)),
    p,
    width = w,
    height = h,
    dpi = 150
  )
}

plot_and_save(
  DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.1", label = TRUE) +
    ggtitle("UMAP res 0.1"),
  "01_umap_res0.1.png"
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.2", label = TRUE) +
    ggtitle("UMAP res 0.2"),
  "02_umap_res0.2.png"
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.5", label = TRUE) +
    ggtitle("UMAP res 0.5"),
  "03_umap_res0.5.png"
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.8", label = TRUE) +
    ggtitle("UMAP res 0.8"),
  "04_umap_res0.8.png"
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "sample", pt.size = 0.3) +
    ggtitle("UMAP by Sample"),
  "05_umap_sample.png",
  w = 14
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "tissue", pt.size = 0.5) +
    ggtitle("UMAP by Tissue"),
  "06_umap_tissue.png"
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "timepoint", pt.size = 0.3) +
    ggtitle("UMAP by Timepoint"),
  "07_umap_timepoint.png",
  w = 14
)

plot_and_save(
  DimPlot(seurat_harmony, group.by = "donor", pt.size = 0.3) +
    ggtitle("UMAP by Donor"),
  "08_umap_donor.png",
  w = 14
)

message("Saving RDS...")
saveRDS(seurat_harmony, paste0("results/", STUDY, "_seurat_processed.rds"))

message("Done.")
