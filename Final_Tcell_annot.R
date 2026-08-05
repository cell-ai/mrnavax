################################################################################
##                                                                            ##
##   T CELL MARKER EXPLORATION — USING THE EXISTING T-CELL-ONLY OBJECT      ##
##                                                                            ##
##   Your cell_type labels here were TRANSFERRED from the main object, so    ##
##   they're the same coarse categories (naive, CD4 em, CD8 eff, ...) — they ##
##   carry no signal one way or the other about Th1/Th2/Th17/Treg/Tscm       ##
##   structure. This script:                                                 ##
##     1. Checks what embedding you've inherited (full-object PCA/UMAP, or   ##
##        something computed on this T-cell object specifically)            ##
##     2. Scores every signature per cell                                   ##
##     3. Visualizes scores on the EXISTING embedding + against the         ##
##        transferred labels — a first look, before committing to a redo    ##
##     4. Gives you a toggle to recompute PCA/UMAP on this object if the     ##
##        first look doesn't show usable structure                          ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(patchwork)
})


# ─────────────────────────────────────────────────────────────────────────────
# 0a. READ IN THE T-CELL OBJECT FROM .h5ad
#     h5ad (AnnData) is the Python/Scanpy container format — converting to
#     Seurat needs one extra hop. Using zellkonverter here (Bioconductor):
#     it manages its own Python/conda env internally via basilisk, so it
#     doesn't need a pre-existing reticulate/scanpy setup to work.
#
#     If you already have tcell_obj loaded as an .rds, skip this block —
#     it only runs when tcell_obj doesn't already exist in the environment.
# ─────────────────────────────────────────────────────────────────────────────
setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")
H5AD_PATH <- "final_seurat_objs/T_scgpt.h5ad"

if (!exists("tcell_obj")) {
  
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("zellkonverter", update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages({
    library(zellkonverter)
    library(SingleCellExperiment)
  })
  
  cat("── Reading h5ad via zellkonverter ──\n")
  sce <- readH5AD(H5AD_PATH)
  
  cat("Assays found in .h5ad (X / layers):", paste(assayNames(sce), collapse = ", "), "\n")
  cat("Embeddings found in .h5ad (obsm):   ", paste(reducedDimNames(sce), collapse = ", "), "\n\n")
  
  # Your object has NO "X" assay — anndata.X was empty and the real data
  # lives in layers "counts" (raw) and "logcounts" (normalized). Point
  # as.Seurat() at those instead:
  tcell_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
  
  # as.Seurat() carries reducedDims(sce) over into tcell_obj@reductions.
  # Your object has several: HARMONY, PCA, UMAP, UMAP_SUB, X_scGPT, X_umap.
  #   - UMAP_SUB is the one to look at FIRST — the "_SUB" naming strongly
  #     suggests it was computed on the T-cell subset itself (as opposed to
  #     inherited from the full multi-lineage object), which is exactly the
  #     kind of embedding that can actually resolve Th1/Th2/Th17/Treg/etc.
  #     structure. UMAP (no suffix) and X_umap are more likely holdovers
  #     from the full-object embedding — good to compare against, not to
  #     lead with.
  #   - HARMONY is a batch-corrected PCA-style space. If you do end up
  #     needing to recompute clusters (Section 4), use HARMONY dims for
  #     FindNeighbors/RunUMAP instead of a fresh RunPCA — it already has
  #     batch effects removed, which a plain re-run of PCA on this subset
  #     would not.
  #   - X_scGPT is a foundation-model (scGPT) per-cell embedding — another
  #     candidate low-dim space if HARMONY/PCA don't resolve well, though
  #     start with UMAP_SUB/HARMONY first.
  reduction_names <- Reductions(tcell_obj)
  cat("Reductions carried over into the Seurat object:", paste(reduction_names, collapse = ", "), "\n\n")
  
  # obs columns (scanpy .obs, including your transferred cell_type) land
  # directly in tcell_obj@meta.data — sanity check the label column name:
  if (!"cell_type" %in% colnames(tcell_obj@meta.data)) {
    cat("[WARN] No 'cell_type' column found in meta.data. Available columns:\n")
    print(colnames(tcell_obj@meta.data))
    cat("Rename whichever column holds your transferred labels to 'cell_type',",
        "or update the group_by calls later in this script.\n\n")
  }
  
  # If NormalizeData was never run in scanpy and X is raw counts, normalize
  # now so AddModuleScore has log-normalized data to work with:
  # tcell_obj <- NormalizeData(tcell_obj, verbose = FALSE)
  
  # Alternatives if zellkonverter is unavailable/too slow to install:
  #   sceasy::convertFormat(H5AD_PATH, from = "anndata", to = "seurat")
  #   schard::h5ad2seurat(H5AD_PATH)   # remotes::install_github("cellgeni/schard")
}


# ─────────────────────────────────────────────────────────────────────────────
# 0.  WHAT EMBEDDING DID WE INHERIT?
#     If this object was made with subset() and nobody reran RunPCA/RunUMAP
#     afterward, the PCA loadings were fit on variance across the WHOLE
#     dataset (B cell / myeloid / T cell separation dominates), and T cells
#     will likely sit as one relatively flat blob on it — even if real
#     Th1/Th2/Th17/Treg structure exists in the underlying expression data.
#     A blob here means "this embedding isn't built to show it," not
#     "the structure doesn't exist." Check before drawing conclusions.
# ─────────────────────────────────────────────────────────────────────────────

cat("── Reductions present on this object ──\n")
print(Reductions(tcell_obj))

cat("\n── Variable features ──\n")
n_vf <- length(VariableFeatures(tcell_obj))
cat("n =", n_vf, "\n")
if (n_vf == 0 && exists("sce")) {
  # as.Seurat() doesn't carry scanpy's .var.highly_variable flags over on
  # its own — recover them directly from the SCE's row metadata instead,
  # so you get the actual scanpy-computed HVG list rather than nothing.
  hvg_col <- intersect(c("highly_variable", "highly_variable_features", "hvg"),
                       colnames(rowData(sce)))
  if (length(hvg_col) > 0) {
    hvg_col   <- hvg_col[1]
    hvg_genes <- rownames(sce)[as.logical(rowData(sce)[[hvg_col]])]
    VariableFeatures(tcell_obj) <- hvg_genes
    cat("Recovered", length(hvg_genes), "HVGs from anndata.var$", hvg_col, "\n\n")
  } else {
    cat("No highly_variable-style column in anndata.var. Columns available:",
        paste(colnames(rowData(sce)), collapse = ", "), "\n")
    cat("(this just means HVGs weren't flagged/exported from scanpy for this object —",
        "not evidence either way about which embedding is lineage-level vs",
        "T-cell-specific; Section 4's fresh-PCA fallback calls",
        "FindVariableFeatures itself regardless, so this isn't a blocker)\n\n")
  }
} else if (n_vf == 0) {
  cat("(0 is expected if tcell_obj came from .rds rather than this session's",
      ".h5ad conversion — `sce` isn't in scope to recover flags from)\n\n")
}

cat("── Transferred cell_type distribution ──\n")
print(table(tcell_obj$cell_type))
cat("\n")

umap_name <- Reductions(tcell_obj)[grepl("umap_sub", Reductions(tcell_obj), ignore.case = TRUE)][1]
if (is.na(umap_name)) {
  cat("No UMAP_SUB found — falling back to a general 'umap' match.\n")
  umap_name <- Reductions(tcell_obj)[grepl("umap", Reductions(tcell_obj), ignore.case = TRUE)][1]
} else {
  cat("Using UMAP_SUB — likely computed on this T-cell subset specifically.\n")
}
if (!is.na(umap_name)) {
  p0 <- DimPlot(tcell_obj, reduction = umap_name, group.by = "cell_type", label = TRUE) +
    ggtitle(paste("Reduction:", umap_name, "— colored by TRANSFERRED cell_type"))
  print(p0)
} else {
  cat("No UMAP reduction found on this object — skipping Step 0 DimPlot.\n")
  cat("(if there's genuinely no reduction at all, you'll need to run RunPCA/RunUMAP",
      "at least once before any of the visualization below is meaningful)\n\n")
}


# ─────────────────────────────────────────────────────────────────────────────
# 1.  HUMAN-ORTHOLOG MARKER SIGNATURES  (unchanged from before)
# ─────────────────────────────────────────────────────────────────────────────

IEL_sig        <- c("CCL5", "NKG7", "ITGAE", "CD8A", "JAML", "GZMB", "GZMA", "CD7", "CD244")
TFH_sig        <- c("BCL6", "CXCR5", "IL21", "PDCD1", "ICOS", "CXCR4")
Th1_sig        <- c("TBX21", "ID2", "IFNG", "IL2RB", "STAT1", "CXCR3", "CCR5")
Th2_sig        <- c("GATA3", "IRF4", "STAT5A", "STAT5B", "IL1RL1", "CXCR4", "IL4")
Th17_sig       <- c("IL17A", "IL17F", "RORC", "RORA", "CCR6", "IL23R", "IL6R", "STAT3")
Treg_sig       <- c("FOXP3", "CTLA4", "TNFRSF18", "IKZF2", "TIGIT", "IZUMO1", "MAF", "IL2RA")
Tscm_sig       <- c("TCF7", "LEF1", "EOMES", "FOXP1", "CERS6", "BCL2", "FOXO1")
Tcmp_sig       <- c("CXCR5", "CCR7", "BCL2", "ID3", "KLF2", "TCF7", "S100A6", "SLAMF6")
Tmem_sig       <- c("IL7R", "CCR7", "BACH2", "CCND2", "EIF4A2")
Tnaive_sig     <- c("SELL", "TCF7", "LEF1", "BACH2", "DAPL1", "KLF2", "CCR7", "S1PR1")
Cytotoxic_sig  <- c("GZMB", "GZMA", "IFNG", "PRF1", "CD8A", "TBX21", "PRDM1")
Cycling_sig    <- c("CCNA2", "CENPF", "CCNB2", "CDKN3", "TPX2", "CDCA3", "PCLAF")
gut_homing_sig <- c("CCR9", "ITGA4", "ITGB7", "ITGB2")

t_cell_signatures <- list(
  IEL = IEL_sig, TFH = TFH_sig, Th1 = Th1_sig, Th2 = Th2_sig, Th17 = Th17_sig,
  Treg = Treg_sig, Tscm = Tscm_sig, Tcmp = Tcmp_sig, Tmem = Tmem_sig,
  Tnaive = Tnaive_sig, Cytotoxic = Cytotoxic_sig, Cycling = Cycling_sig,
  gut_homing = gut_homing_sig
)

missing_genes <- lapply(t_cell_signatures, function(g) setdiff(g, rownames(tcell_obj)))
missing_genes <- missing_genes[lengths(missing_genes) > 0]
if (length(missing_genes) > 0) {
  cat("── Genes NOT found in object (dropped by AddModuleScore) ──\n")
  print(missing_genes)
  cat("\n")
}


# ─────────────────────────────────────────────────────────────────────────────
# 2.  SCORE EVERY SIGNATURE PER CELL
# ─────────────────────────────────────────────────────────────────────────────

tcell_obj <- AddModuleScore(tcell_obj, features = t_cell_signatures, name = "sig_")
score_cols <- paste0("sig_", seq_along(t_cell_signatures))
sig_names  <- names(t_cell_signatures)
colnames(tcell_obj@meta.data)[match(score_cols, colnames(tcell_obj@meta.data))] <- sig_names


# ─────────────────────────────────────────────────────────────────────────────
# 3.  FIRST LOOK — scores on the existing embedding + against transferred labels
#     This is the "let's just look" step you asked for, BEFORE deciding
#     whether a recompute is needed.
# ─────────────────────────────────────────────────────────────────────────────

if (!is.na(umap_name)) {
  p1 <- FeaturePlot(tcell_obj, features = sig_names, reduction = umap_name,
                    ncol = 4, min.cutoff = "q05", max.cutoff = "q95")
  print(p1)
  ggsave("tcell_signature_umaps_EXISTING_embedding.pdf", p1, width = 16, height = 12)
}

# Tabular check independent of embedding geometry: does any TRANSFERRED
# cell_type already skew toward a particular signature?
score_by_label <- tcell_obj@meta.data %>%
  group_by(RNA_snn_res.0.5) %>%
  summarise(across(all_of(sig_names), mean), n_cells = n(), .groups = "drop")

cat("── Mean signature score per TRANSFERRED cell_type ──\n")
print(score_by_label, width = Inf, n = Inf)

mat <- as.matrix(score_by_label[, sig_names])
rownames(mat) <- score_by_label$RNA_snn_res.0.5

pdf("tcell_signature_heatmap_by_transferred_label.pdf", width = 8, height = 6)
Heatmap(t(scale(mat)), name = "z-score",
        column_title = "Mean signature score per transferred cell_type")
dev.off()


# ─────────────────────────────────────────────────────────────────────────────
# UNBIASED MARKER GENES PER CLUSTER (resolution 0.8)
#   - Confirms/refines the signature-based labels above
#   - Reveals what's actually driving ambiguous clusters (5, 22) that the
#     13 pre-defined signatures didn't cleanly resolve
#   - Excludes near-singleton clusters (<50 cells) from the test, since
#     those are noise/doublets, not real groups worth testing
# ─────────────────────────────────────────────────────────────────────────────

Idents(tcell_obj) <- "RNA_snn_res.0.5"   # match the column you clustered on
DefaultAssay(tcell_obj) <- "originalexp"

cluster_sizes   <- table(Idents(tcell_obj))
tiny_clusters   <- names(cluster_sizes)[cluster_sizes < 200]
cat("Excluding likely-noise clusters (<50 cells):", paste(tiny_clusters, collapse = ", "), "\n")

tcell_obj_clean <- subset(tcell_obj, idents = setdiff(levels(Idents(tcell_obj)), tiny_clusters))

markers_res0.5 <- FindAllMarkers(
  tcell_obj_clean,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

top_markers <- markers_res0.5 %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 30) %>%
  ungroup()

write.csv(markers_res0.8, "tcell_findallmarkers_res0.8_full.csv", row.names = FALSE)
write.csv(top_markers,    "tcell_findallmarkers_res0.8_top15.csv", row.names = FALSE)
print(top_markers, n = Inf)

# ─────────────────────────────────────────────────────────────────────────────
# LABEL CLUSTERS FOR DimPlot (resolution 0.5, cluster 14/B-cells already removed)
#   Clusters 0 and 2 merged — same core markers (TMIGD2, NOSIP, BACH2, SATB1),
#   only difference was RPS4Y1 (sex-linked, not distinct biology).
#   Cluster 8 kept SEPARATE and labeled as technical — pure IEG/stress program
#   with no lineage marker; merging it into a real cell type would misrepresent
#   what those cells actually are.
#   Cluster 9 kept separate too, flagged low-confidence pending the n_cells
#   check — merge into "Naive T cells" yourself if it turns out to be tiny/
#   sex-linked like 2 was.
# ─────────────────────────────────────────────────────────────────────────────

CLUSTER_COL <- "RNA_snn_res.0.5"   # change to whatever column you clustered on

cluster_labels <- c(
  "0"  = "Naive T cells",
  "1"  = "Memory T cells",
  "2"  = "Naive T cells",                          # merged with 0
  "3"  = "CD8+ T cells (NKG2+ cytotoxic)",
  "4"  = "Activated CD4 T helper cells",
  "5"  = "Regulatory T cells (Treg)",
  "6"  = "T follicular helper (Tfh)",
  "7"  = "Activated/effector T cells",
  "8"  = "Stress-response signature (technical)",
  "9"  = "Low-confidence / sex-linked cluster",
  "10" = "IFN-stimulated T cells",
  "11" = "Cytotoxic/effector-memory T cells (GZMK+)",
  "12" = "Tr1 cells (IL10+ regulatory)",
  "13" = "IFN-stimulated CD8+ T cells"
)

new_labels <- cluster_labels[as.character(tcell_obj[[CLUSTER_COL]][, 1])]

tcell_obj$cell_type_annotated <- factor(
  unname(new_labels),
  levels = unique(cluster_labels)
)

table(tcell_obj$cell_type_annotated)

# Visualize on UMAP_SUB (or whichever reduction you've been using)
sub_umap <- Reductions(tcell_obj)[grepl("umap_sub", Reductions(tcell_obj), ignore.case = TRUE)][1]

DimPlot(tcell_obj, reduction = sub_umap, group.by = "cell_type_annotated",
        label = TRUE, repel = TRUE) +
  ggtitle("T cell subsets (resolution 0.5, annotated)")


# ─────────────────────────────────────────────────────────────────────────────
# EXCLUDE SEX-LINKED GENES, THEN RE-RUN FindAllMarkers
#   Removes them from the TEST, not from clustering — clusters 2/9 will
#   still exist as-is, but this reveals whether their real (non-sex-gene)
#   markers look distinct/biological or just echo another cluster, which
#   is the actual evidence needed to decide whether to merge them.
# ─────────────────────────────────────────────────────────────────────────────

sex_linked_genes <- c(
  # Y chromosome (male-specific)
  "RPS4Y1", "RPS4Y2", "DDX3Y", "KDM5D", "UTY", "USP9Y", "ZFY",
  "EIF1AY", "TXLNGY", "NLGN4Y", "TMSB4Y", "TTTY14", "TTTY15",
  # X-inactivation-associated (expressed almost exclusively from Xi)
  "XIST", "TSIX", "FTX", "JPX"
)

present_sex_genes <- intersect(sex_linked_genes, rownames(tcell_obj))
cat("Sex-linked genes found in object (excluded from test):",
    paste(present_sex_genes, collapse = ", "), "\n")
cat("Not present (nothing to exclude for these):",
    paste(setdiff(sex_linked_genes, present_sex_genes), collapse = ", "), "\n\n")

features_to_test <- setdiff(rownames(tcell_obj_clean), present_sex_genes)

Idents(tcell_obj_clean) <- "RNA_snn_res.0.5"   # match whatever you clustered on

markers_no_sex <- FindAllMarkers(
  tcell_obj_clean,
  features        = features_to_test,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

top_markers_no_sex <- markers_no_sex %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

write.csv(markers_no_sex,     "tcell_findallmarkers_no_sexgenes_full.csv", row.names = FALSE)
write.csv(top_markers_no_sex, "tcell_findallmarkers_no_sexgenes_top15.csv", row.names = FALSE)
print(top_markers_no_sex, n = Inf)

CLUSTER_COL <- "RNA_snn_res.0.5"

# Drop the two negligible clusters before reassignment so they can't
# contaminate the centroids or the reassignment target set
tcell_obj_final <- subset(tcell_obj, subset = !(RNA_snn_res.0.5 %in% c("14", "15")))

technical_clusters <- c("8", "9")
real_clusters <- setdiff(levels(tcell_obj_final[[CLUSTER_COL]][, 1]),
                         c(technical_clusters, "14", "15"))

real_centroids <- tcell_obj_final@meta.data %>%
  filter(.data[[CLUSTER_COL]] %in% real_clusters) %>%
  group_by(.data[[CLUSTER_COL]]) %>%
  summarise(across(all_of(sig_names), mean), .groups = "drop")

centroid_mat <- as.matrix(real_centroids[, sig_names])
rownames(centroid_mat) <- as.character(real_centroids[[CLUSTER_COL]])

technical_mask   <- tcell_obj_final[[CLUSTER_COL]][, 1] %in% technical_clusters
technical_scores <- as.matrix(tcell_obj_final@meta.data[technical_mask, sig_names])

nearest_real_cluster <- apply(technical_scores, 1, function(cell_vec) {
  cors <- apply(centroid_mat, 1, function(centroid_vec) cor(cell_vec, centroid_vec))
  names(cors)[which.max(cors)]
})

final_cluster <- as.character(tcell_obj_final[[CLUSTER_COL]][, 1])
names(final_cluster) <- colnames(tcell_obj_final)
final_cluster[technical_mask] <- nearest_real_cluster

tcell_obj_final$cluster_reassigned <- final_cluster
table(nearest_real_cluster)   # <- paste this so we can sanity-check the redistribution

cluster_labels <- c(
  "0"  = "Naive T cells",
  "1"  = "Memory T cells",
  "2"  = "Naive T cells",
  "3"  = "CD8+ T cells (NKG2+ cytotoxic)",
  "4"  = "Activated CD4 T helper cells",
  "5"  = "Regulatory T cells (Treg)",
  "6"  = "T follicular helper (Tfh)",
  "7"  = "Activated/effector T cells",
  "10" = "IFN-stimulated T cells",
  "11" = "Cytotoxic/effector-memory T cells (GZMK+)",
  "12" = "Tr1 cells (IL10+ regulatory)",
  "13" = "IFN-stimulated CD8+ T cells"
)

tcell_obj_final$cell_type_annotated <- factor(
  unname(cluster_labels[tcell_obj_final$cluster_reassigned]),
  levels = unique(cluster_labels)
)

table(tcell_obj_final$cell_type_annotated)

DimPlot(tcell_obj_final, reduction = sub_umap, group.by = "cell_type_annotated",
        label = TRUE, repel = TRUE) + ggtitle("T cell subsets — final annotation")

# ─────────────────────────────────────────────────────────────────────────────
# ISOLATE THE Th (CD4 helper) CLUSTER AND SUBCLUSTER FOR Th1/Th2/Th17
# ─────────────────────────────────────────────────────────────────────────────
VlnPlot(tcell_obj_final, features = c("CD4", "CD8A", "CD8B", "CD3E"),
        idents = "4", group.by = "cluster_reassigned")

saveRDS(tcell_obj_final, "final_seurat_objs/Tcell_w_final_annot.RDS")

################################################################################
##                                                                            ##
##   Th SUBCLUSTERING: ISOLATE, RECLUSTER, SCORE, ANNOTATE                  ##
##                                                                            ##
##   Starting point: tcell_obj_final, with a working whole-object cluster    ##
##   annotation already in cluster_reassigned (from the earlier T-cell       ##
##   pipeline: technical clusters 8/9 reassigned by nearest signature        ##
##   centroid, clusters 14/15 dropped for being too small to trust).         ##
##                                                                            ##
##   This script:                                                            ##
##     1. Isolates the CD4 Th cluster (cluster "4" in cluster_reassigned)    ##
##     2. Re-embeds/re-clusters on JUST this subset (res 0.2)                ##
##     3. Scores Th1/Th2/Th17 signatures with subset-specific names          ##
##        (avoids column collisions with the whole-object AddModuleScore     ##
##        run from earlier)                                                  ##
##     4. Runs FindAllMarkers + a canonical-marker DotPlot, and pulls the    ##
##        DotPlot's underlying numbers into a clean table                   ##
##     5. Reassigns technical Th subclusters to their nearest real neighbor  ##
##        (same logic as the whole-object technical-cluster handling)       ##
##     6. Applies final labels based on the marker/DotPlot evidence         ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})


# ─────────────────────────────────────────────────────────────────────────────
# 1. ISOLATE THE Th (CD4 HELPER) CLUSTER FROM THE REST OF THE T-CELL OBJECT
# ─────────────────────────────────────────────────────────────────────────────

# Confirm lineage before isolating — this is the check that told us cluster 4
# is genuinely CD4+ (broad CD4/CD3E, near-zero CD8A/CD8B) rather than assuming it
VlnPlot(tcell_obj_final, features = c("CD4", "CD8A", "CD8B", "CD3E"),
        idents = "4", group.by = "cluster_reassigned")

th_obj <- subset(tcell_obj_final, subset = cluster_reassigned == "4")
cat("Th cluster isolated:", ncol(th_obj), "cells\n\n")


# ─────────────────────────────────────────────────────────────────────────────
# 2. RE-EMBED / RE-CLUSTER ON THE Th SUBSET SPECIFICALLY
#    Removing naive/memory/Treg/Tfh/CD8 identity (the dominant variance axes
#    in the whole T-cell object) lets the much weaker Th1/Th2/Th17 programs
#    actually drive the PCA/clustering here.
# ─────────────────────────────────────────────────────────────────────────────

harmony_name <- Reductions(th_obj)[grepl("harmony", Reductions(th_obj), ignore.case = TRUE)][1]

if (!is.na(harmony_name)) {
  n_dims <- ncol(Embeddings(th_obj, harmony_name))
  th_obj <- FindNeighbors(th_obj, reduction = harmony_name, dims = 1:min(30, n_dims), verbose = FALSE)
  umap_reduction <- harmony_name
} else {
  th_obj <- NormalizeData(th_obj, verbose = FALSE)
  th_obj <- FindVariableFeatures(th_obj, nfeatures = 2000, verbose = FALSE)
  th_obj <- ScaleData(th_obj, verbose = FALSE)
  th_obj <- RunPCA(th_obj, npcs = 30, verbose = FALSE)
  th_obj <- FindNeighbors(th_obj, dims = 1:30, verbose = FALSE)
  umap_reduction <- "pca"
}

RESOLUTION <- 0.2   # settled on this after comparing 0.6 vs 0.2 — 0.2 gave the
# cleaner, more reproducible split (RORC/CCR6/KLRB1 all
# converged on one cluster at this resolution)

th_obj <- FindClusters(th_obj, resolution = RESOLUTION, verbose = FALSE)
th_obj <- RunUMAP(th_obj, reduction = umap_reduction, dims = 1:30, verbose = FALSE)

cat("Th subclusters found:", n_distinct(th_obj$seurat_clusters), "\n")
print(table(th_obj$seurat_clusters))
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# 3. SCORE Th1/Th2/Th17 ON THE SUBSET (distinct column names — the object
#    already carries whole-object Th1/Th2/Th17 columns from the earlier
#    13-signature AddModuleScore run; reusing those names would collide)
# ─────────────────────────────────────────────────────────────────────────────

th_sigs <- list(Th1 = Th1_sig, Th2 = Th2_sig, Th17 = Th17_sig)   # defined earlier in the pipeline

# Guard against re-running this section on an object that already has these
# columns (from a previous run of this script, or an earlier AddModuleScore
# call) — AddModuleScore's rename step appends rather than overwrites, which
# causes a duplicate-column error on group_by/summarise if not cleared first.
th_obj@meta.data[, intersect(c("Th1_sub", "Th2_sub", "Th17_sub"),
                             colnames(th_obj@meta.data))] <- NULL

th_obj <- AddModuleScore(th_obj, features = th_sigs, name = "th_sub_")
colnames(th_obj@meta.data)[match(paste0("th_sub_", 1:3), colnames(th_obj@meta.data))] <-
  c("Th1_sub", "Th2_sub", "Th17_sub")

score_by_subcluster <- th_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(c(Th1_sub, Th2_sub, Th17_sub), mean), n_cells = n(), .groups = "drop")
cat("── Composite signature scores per subcluster (diluted by cytokine dropout —\n",
    "   treat as a rough cross-check only, NOT primary evidence; see Section 4) ──\n")
print(as.data.frame(score_by_subcluster))
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# 4. UNBIASED MARKERS + CANONICAL MASTER-TF/SURFACE-MARKER CHECK
#    This is the primary evidence — composite scores get diluted by cytokine
#    dropout and coincidental single-gene overlaps; individual canonical
#    genes (TBX21, GATA3, RORC + their surface markers) are more reliable.
# ─────────────────────────────────────────────────────────────────────────────

Idents(th_obj) <- "seurat_clusters"
th_markers <- FindAllMarkers(th_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
th_top <- th_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15) %>% ungroup()
print(th_top, n = Inf)

# Canonical panel — drop clusters already identified as technical once you
# know which those are for your run (see Section 5); shown here unfiltered
# the first time through so you can spot the technical ones from the markers.
th_canonical_genes <- c(
  "TBX21", "CXCR3", "STAT1",              # Th1
  "GATA3", "IL4R", "CCR4",                # Th2
  "RORC", "RORA", "CCR6", "KLRB1",        # Th17
  "CXCR5", "PDCD1", "TOX2", "IL21"        # Tfh (checks whether a "Th2" cluster is actually Tfh-skewed)
)

p_dot <- DotPlot(th_obj, features = th_canonical_genes, group.by = "seurat_clusters")

dotplot_table <- p_dot$data %>%
  dplyr::select(cluster = id, gene = features.plot,
                avg_exp_scaled = avg.exp.scaled, avg_exp = avg.exp,
                pct_expressed = pct.exp)

scaled_wide <- dotplot_table %>%
  dplyr::select(cluster, gene, avg_exp_scaled) %>%
  pivot_wider(names_from = gene, values_from = avg_exp_scaled)

cat("── Average expression (scaled) per cluster, canonical Th markers ──\n")
print(as.data.frame(scaled_wide))
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# 5. FINAL CHECKS before committing to labels
# ─────────────────────────────────────────────────────────────────────────────

# (a) Is the "Th2/Tfh-skewed" cluster one mixed population, or two
#     intermixed cell types this resolution didn't separate?
#     -> replace "4" with whichever cluster ID showed high GATA3+TOX2
#        together in your scaled_wide table
FeatureScatter(th_obj, feature1 = "TOX2", feature2 = "CXCR5",
               cells = WhichCells(th_obj, idents = "4")) +
  ggtitle("Candidate Th2/Tfh-skewed cluster: TOX2 vs CXCR5, per cell")

# (b) Does the Th17 call hold up on a broader (if sparser) cytokine panel?
#     -> replace "2" with whichever cluster ID showed RORC+RORA+CCR6+KLRB1
VlnPlot(th_obj, features = c("IL17A", "IL17F", "IL22", "AHR"),
        idents = "2", ncol = 4)


# ─────────────────────────────────────────────────────────────────────────────
# 6. REASSIGN TECHNICAL SUBCLUSTERS, THEN APPLY FINAL LABELS
#    Same logic as the whole-object technical clusters: reassign per cell
#    to the nearest REAL Th subcluster by signature-score correlation,
#    rather than bulk-merging a heterogeneous technical cluster into one
#    real cell type.
# ─────────────────────────────────────────────────────────────────────────────

# Update these once your run confirms which subcluster IDs are technical vs
# real — based on the walkthrough that produced this script, technical
# clusters showed long-transcript/no-lineage-marker signatures (NEAT1, ANK3,
# SYNE2, MACF1...) or pure dissociation-stress genes (NR4A1, EGR1, FOS/JUN...)
th_technical_clusters <- c("3", "5")
th_real_clusters <- setdiff(levels(th_obj$seurat_clusters), th_technical_clusters)

real_centroids <- th_obj@meta.data %>%
  filter(seurat_clusters %in% th_real_clusters) %>%
  group_by(seurat_clusters) %>%
  summarise(across(c(Th1_sub, Th2_sub, Th17_sub), mean), .groups = "drop")

centroid_mat <- as.matrix(real_centroids[, c("Th1_sub", "Th2_sub", "Th17_sub")])
rownames(centroid_mat) <- as.character(real_centroids$seurat_clusters)

technical_mask   <- th_obj$seurat_clusters %in% th_technical_clusters
technical_scores <- as.matrix(th_obj@meta.data[technical_mask, c("Th1_sub", "Th2_sub", "Th17_sub")])

nearest_real_cluster <- apply(technical_scores, 1, function(cell_vec) {
  cors <- apply(centroid_mat, 1, function(centroid_vec) cor(cell_vec, centroid_vec))
  names(cors)[which.max(cors)]
})

th_final_cluster <- as.character(th_obj$seurat_clusters)
names(th_final_cluster) <- colnames(th_obj)
th_final_cluster[technical_mask] <- nearest_real_cluster
th_obj$th_cluster_reassigned <- th_final_cluster

cat("── Technical Th subcluster reassignment ──\n")
print(table(nearest_real_cluster))
cat("\n")

# Labels based on the canonical-marker evidence (Section 4):
#   - RORC/RORA/CCR6/KLRB1 all converge, Th1/Th2-low  -> Th17
#   - GATA3/IL4R/CCR4/TOX2 high, but CXCR5 low         -> Th2, partial Tfh features
#   - No cluster showed TBX21 + downstream targets together -> no discrete Th1 found
th_cluster_labels <- c(
  "0" = "Quiescent/resting Th cells",
  "1" = "Activated Th cells (no clear lineage skew)",
  "2" = "Th17 cells",
  "4" = "Th2 cells (GATA3-high, partial Tfh-like features)"
)

th_obj$th_annotated <- factor(
  unname(th_cluster_labels[th_obj$th_cluster_reassigned]),
  levels = unique(th_cluster_labels)
)

table(th_obj$th_annotated)

DimPlot(th_obj, group.by = "th_annotated", label = TRUE, repel = TRUE) +
  ggtitle("Th subsets — final annotation")

FeaturePlot(th_obj, features = c("Th1_sub", "Th2_sub", "Th17_sub"),
            ncol = 3, min.cutoff = "q05", max.cutoff = "q95") +
  plot_annotation(title = "Th1/Th2/Th17 composite scores on Th subset UMAP")
