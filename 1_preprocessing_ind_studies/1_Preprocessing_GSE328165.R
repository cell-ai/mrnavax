#!/usr/bin/env Rscript
# =============================================================
# 1.Processing.R
# Owner: Izabela Mamede
# Date:  2026-05-17
# Altered by Marcela
# =============================================================

setwd("/media/csbl/sandbox-SSD-3/mra_vaccine")

STUDY <- "GSE328165"

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(harmony)
})

# ── Parallel settings ─────────────────────────────────────────

# ── Parallel settings ─────────────────────────────────────────
N_THREADS <- 10
options(future.globals.maxSize = 200 * 1024^3)   # 200 GB — covers large objects
future::plan("multicore", workers = N_THREADS)

message("Running with ", N_THREADS, " threads")

dir.create("results/qc_plots",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/umap_plots", recursive = TRUE, showWarnings = FALSE)

# ── Metadata helpers ──────────────────────────────────────────

# Generic: return first regex match or a fallback string
extract_or_unknown <- function(x, pattern, fallback = "unknown") {
  m <- regmatches(x, regexpr(pattern, x, perl = TRUE))
  if (length(m) == 0 || is.na(m) || nchar(m) == 0) fallback else m
}

# Tissue: check both _PB_ / _PBMC_ and _LN_ / _lymph_ variants
parse_tissue <- function(lib) {
  dplyr::case_when(
    grepl("_PB_|_PBMC_",          lib, ignore.case = TRUE) ~ "PBMC",
    grepl("_LN_|_lymph",          lib, ignore.case = TRUE) ~ "lymph_node",
    grepl("blood",                 lib, ignore.case = TRUE) ~ "PBMC",
    grepl("lymph|node|tonsil",     lib, ignore.case = TRUE) ~ "lymph_node",
    TRUE ~ "unknown"
  )
}

# Timepoint: matches d0, d7, d14, day0, day7, wk1, w1, etc.
parse_timepoint <- function(lib) {
  tp <- extract_or_unknown(lib, "(?i)(day|d|wk|w)\\d+", fallback = NA_character_)
  if (is.na(tp)) {
    # last-resort: bare number preceded by nothing letter-like
    tp <- extract_or_unknown(lib, "(?<=[-_])\\d{1,3}(?=[-_]|$)", fallback = "unknown")
  }
  tolower(tp)
}

# Donor: ELAB-WU### or WU### or donor### or Pt### or S###, etc.
parse_donor <- function(lib) {
  patterns <- c(
    "ELAB-(WU\\d+)"  = "\\1",   # ELAB-WU01  -> WU01
    "(WU\\d+)"       = "\\1",   # WU01       -> WU01
    "(donor\\d+)"    = "\\1",   # donor3     -> donor3
    "(Pt\\d+)"       = "\\1",   # Pt5        -> Pt5
    "(\\bS\\d{2,})"  = "\\1"    # S07        -> S07
  )
  for (pat in names(patterns)) {
    m <- regmatches(lib, regexpr(pat, lib, perl = TRUE))
    if (length(m) > 0 && nchar(m) > 0)
      return(gsub(pat, patterns[[pat]], m, perl = TRUE))
  }
  "unknown"
}

# ── Core splitter ─────────────────────────────────────────────

split_aggr_into_samples <- function(matrix_path, features_path,
                                    barcodes_path, aggregation_csv) {

  message("Reading matrix...")
  mat      <- as(readMM(matrix_path), "CsparseMatrix")
  features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  agg      <- read.csv(aggregation_csv, stringsAsFactors = FALSE)

  # ── Validate aggregation CSV columns ────────────────────────
  required_cols <- c("sample_id", "library_id")
  missing_cols  <- setdiff(required_cols, colnames(agg))
  if (length(missing_cols) > 0) {
    stop("aggregation_csv is missing required column(s): ",
         paste(missing_cols, collapse = ", "),
         "\n  Found: ", paste(colnames(agg), collapse = ", "))
  }

  message("Aggregation CSV contains ", nrow(agg), " samples:")
  print(agg[, intersect(c("sample_id", "library_id", "molecule_h5"), colnames(agg))])

  # ── Gene names ───────────────────────────────────────────────
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

  # ── Barcode → sample index ───────────────────────────────────
  # CellRanger aggr appends "-N" (1-based) at the end of each barcode.
  # Robustly extract that suffix; warn clearly if it is missing.
  barcode_sample_idx <- suppressWarnings(
    as.integer(sub(".*-(\\d+)$", "\\1", barcodes$V1))
  )

  n_na <- sum(is.na(barcode_sample_idx))
  if (n_na > 0) {
    # Show a few examples to help the user debug their CSV / barcodes
    bad_examples <- head(barcodes$V1[is.na(barcode_sample_idx)], 5)
    stop(n_na, " barcodes do not end with a numeric sample suffix (e.g. '-1', '-2').\n",
         "  Examples: ", paste(bad_examples, collapse = ", "), "\n",
         "  Check that the barcodes.tsv.gz comes from a CellRanger aggr run.")
  }

  observed_indices <- sort(unique(barcode_sample_idx))
  if (!all(seq_len(nrow(agg)) %in% observed_indices)) {
    message("  WARNING: some aggregation CSV rows have no barcodes in the matrix.")
    message("  Expected indices: ", paste(seq_len(nrow(agg)), collapse = ", "))
    message("  Observed indices: ", paste(observed_indices, collapse = ", "))
  }

  # ── Build one Seurat object per sample ───────────────────────
  sample_objects <- list()

  for (idx in seq_len(nrow(agg))) {

    sample_id <- agg$sample_id[idx]
    lib       <- agg$library_id[idx]
    cell_mask <- barcode_sample_idx == idx

    if (sum(cell_mask) == 0) {
      message("  [SKIP] ", sample_id,
              " — no cells found for barcode suffix -", idx,
              " (library_id: ", lib, ")")
      next
    }

    sub_mat <- mat[, cell_mask, drop = FALSE]

    obj <- CreateSeuratObject(
      counts       = sub_mat,
      project      = sample_id,
      min.cells    = 3,
      min.features = 200
    )

    obj$sample    <- sample_id
    obj$gse       <- sub("_aggregation.csv$", "", basename(aggregation_csv))
    obj$tissue    <- parse_tissue(lib)
    obj$timepoint <- parse_timepoint(lib)
    obj$donor     <- parse_donor(lib)

    obj <- RenameCells(obj, add.cell.id = sample_id)

    message(sprintf(
      "  [OK] %-30s | tissue: %-12s | timepoint: %-6s | donor: %-8s | cells: %d",
      sample_id,
      unique(obj$tissue),
      unique(obj$timepoint),
      unique(obj$donor),
      ncol(obj)
    ))

    # Warn if metadata looks incomplete
    if (unique(obj$tissue)    == "unknown") message("    ^ tissue not parsed from: ", lib)
    if (unique(obj$timepoint) == "unknown") message("    ^ timepoint not parsed from: ", lib)
    if (unique(obj$donor)     == "unknown") message("    ^ donor not parsed from: ", lib)

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

# ── Merge ─────────────────────────────────────────────────────

if (length(seurat_list) == 0) {
  stop("No samples were loaded. Check aggregation CSV and barcode suffixes.")
}

message("Merging samples...")

if (length(seurat_list) == 1) {
  seurat_merged <- seurat_list[[1]]
} else {
  seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1], project = STUDY)
}

if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])
}

# ── QC ────────────────────────────────────────────────────────

message("Computing QC metrics...")
seurat_merged[["percent_mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")
seurat_merged$log10GenesPerUMI <- log10(seurat_merged$nFeature_RNA) /
                                  log10(seurat_merged$nCount_RNA)

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
  width = 14, height = 12, dpi = 150
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

# ── Normalise / Dimensionality reduction ──────────────────────

DefaultAssay(seurat_qc) <- "RNA"

message("Normalizing...")
#seurat_qc <- NormalizeData(seurat_qc, verbose = FALSE)
#seurat_qc <- FindVariableFeatures(seurat_qc, nfeatures = 2000, verbose = FALSE)
#seurat_qc <- ScaleData(seurat_qc, verbose = FALSE)
#seurat_qc <- RunPCA(seurat_qc, npcs = 50, verbose = FALSE)

#message("Running Harmony...")
#seurat_harmony <- RunHarmony(
#  seurat_qc,
#  group.by.vars = "sample",
#  verbose = FALSE
#)

seurat_qc <- NormalizeData(seurat_qc, verbose = FALSE)
seurat_qc <- FindVariableFeatures(seurat_qc, nfeatures = 2000, verbose = FALSE)
seurat_qc <- ScaleData(seurat_qc, verbose = FALSE)
seurat_qc <- RunPCA(seurat_qc, npcs = 50, verbose = FALSE)

# Harmony, UMAP, clustering can safely use multiple cores
future::plan("multicore", workers = N_THREADS)

message("Running Harmony...")
seurat_harmony <- RunHarmony(seurat_qc, group.by.vars = "sample", verbose = FALSE)

####

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

# ── Clustering ────────────────────────────────────────────────

message("Clustering...")
for (res in c(0.1, 0.2, 0.5, 0.8)) {
  seurat_harmony <- FindClusters(seurat_harmony, resolution = res, verbose = FALSE)
}

# ── UMAP plots ────────────────────────────────────────────────

message("Saving UMAP plots...")

plot_and_save <- function(p, filename, w = 10, h = 8) {
  ggsave(
    file.path("results/umap_plots", paste0(STUDY, "_", filename)),
    p, width = w, height = h, dpi = 150
  )
}

plot_and_save(DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.1", label = TRUE) + ggtitle("UMAP res 0.1"), "01_umap_res0.1.png")
plot_and_save(DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.2", label = TRUE) + ggtitle("UMAP res 0.2"), "02_umap_res0.2.png")
plot_and_save(DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.5", label = TRUE) + ggtitle("UMAP res 0.5"), "03_umap_res0.5.png")
plot_and_save(DimPlot(seurat_harmony, group.by = "RNA_snn_res.0.8", label = TRUE) + ggtitle("UMAP res 0.8"), "04_umap_res0.8.png")
plot_and_save(DimPlot(seurat_harmony, group.by = "sample",    pt.size = 0.3) + ggtitle("UMAP by Sample"),    "05_umap_sample.png",    w = 14)
plot_and_save(DimPlot(seurat_harmony, group.by = "tissue",    pt.size = 0.5) + ggtitle("UMAP by Tissue"),    "06_umap_tissue.png")
plot_and_save(DimPlot(seurat_harmony, group.by = "timepoint", pt.size = 0.3) + ggtitle("UMAP by Timepoint"), "07_umap_timepoint.png", w = 14)
plot_and_save(DimPlot(seurat_harmony, group.by = "donor",     pt.size = 0.3) + ggtitle("UMAP by Donor"),     "08_umap_donor.png",     w = 14)

message("Saving RDS...")
saveRDS(seurat_harmony, paste0("results/", STUDY, "_seurat_processed.rds"))

message("Done.")
