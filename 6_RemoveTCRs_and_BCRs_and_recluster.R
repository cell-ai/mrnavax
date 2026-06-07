setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(DESeq2)
library(volcano3D)
library(Seurat)
library(plotly)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(patchwork)

# ============================================================
# Load and add metadata
# ============================================================

seurat_obj <- readRDS("merged_after_reclustering_w_res_dims20.rds")

sra1 <- read.csv("SraRunTable (4).csv", check.names = FALSE)
sra2 <- read.csv("SraRunTable (5).csv", check.names = FALSE)
sra  <- bind_rows(sra1, sra2)

library_df <- sra %>%
  filter(`Assay Type` == "RNA-Seq") %>%
  transmute(
    sample         = `Sample Name`,
    library_name   = `Library Name`,
    library_run    = Run,
    library_biorep = Biological_replicate,
    sra_study      = `SRA Study`
  ) %>%
  distinct()

meta  <- seurat_obj[[]]
meta2 <- meta %>%
  tibble::rownames_to_column("cell") %>%
  left_join(library_df, by = "sample") %>%
  tibble::column_to_rownames("cell")

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = meta2[, c("library_name", "library_run", "library_biorep", "sra_study")]
)

meta <- seurat_obj[[]] %>%
  mutate(
    patient_id = str_extract(library_name, "WU397-\\d{3}|WU368-\\d{2}"),
    vaccine = case_when(
      str_starts(library_name, "ELAB-WU368") ~ "covidmRNA",
      patient_id %in% c("WU397-006", "WU397-017", "WU397-028", "WU397-029") ~ "mRNA-1010",
      patient_id %in% c("WU397-005", "WU397-009", "WU397-022") ~ "Fluarix",
      TRUE ~ NA_character_
    )
  )

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = meta[, c("patient_id", "vaccine")]
)

cat("Vaccine counts:\n")
print(table(seurat_obj$vaccine, useNA = "ifany"))

# ============================================================
# Compositional plots
# ============================================================

plot_cluster_composition <- function(meta, variable, cluster_col = "RNA_snn_res.0.1") {
  cluster_meta <- meta %>%
    count(
      cluster  = .data[[cluster_col]],
      variable = .data[[variable]],
      name     = "n_cells"
    )
  ggplot(cluster_meta, aes(x = cluster, y = n_cells, fill = variable)) +
    geom_col(position = "fill") +
    theme_bw() +
    labs(x = "Cluster", y = "Proportion of cells", fill = variable,
         title = paste("Cluster composition by", variable))
}

meta <- seurat_obj[[]]
plot_cluster_composition(meta, "study")
plot_cluster_composition(meta, "sex")
plot_cluster_composition(meta, "vaccine")
plot_cluster_composition(meta, "patient_id")
plot_cluster_composition(meta, "timepoint")

# ============================================================
# Step 1 — Remove ALL TCR and BCR genes in one clean pass
# ============================================================

all_genes <- rownames(seurat_obj)

tcr_genes <- grep("^TR[ABGD][VJC]", all_genes, value = TRUE)
bcr_all   <- grep("^IGH|^IGL|^IGK", all_genes, value = TRUE)

# Keep isotype genes — biologically meaningful for vaccine response
keep_isotype <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4",
                  "IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE")

bcr_remove   <- bcr_all[!bcr_all %in% keep_isotype]
genes_remove <- union(tcr_genes, bcr_remove)
genes_keep   <- all_genes[!all_genes %in% genes_remove]

cat(sprintf("Genes before: %d | after: %d | removed: %d (TCR: %d, BCR: %d)\n",
    length(all_genes), length(genes_keep), length(genes_remove),
    length(tcr_genes), length(bcr_remove)))
cat("Isotype genes kept:", paste(keep_isotype[keep_isotype %in% all_genes], collapse = ", "), "\n")

# ============================================================
# Step 2 — Create filtered object from seurat_obj
# ============================================================

seurat_filt <- seurat_obj[genes_keep, ]
cat("Genes in filtered object:", nrow(seurat_filt), "\n")

# Verify nothing slipped through
cat("TCR remaining:", length(grep("^TR[ABGD][VJC]", rownames(seurat_filt))), "\n")
cat("BCR remaining:", length(grep("^IGH|^IGL|^IGK",  rownames(seurat_filt))), "\n")

# ============================================================
# Step 3 — Normalize, variable features, scale, PCA, cluster
# ============================================================

seurat_filt <- NormalizeData(seurat_filt)
seurat_filt <- FindVariableFeatures(seurat_filt, nfeatures = 3000)

var_features <- VariableFeatures(seurat_filt)
cat("TCR in var features:", length(grep("^TR[ABGD][VJC]", var_features)), "\n")
cat("BCR in var features:", length(grep("^IGH|^IGL|^IGK",  var_features)), "\n")

seurat_filt <- ScaleData(seurat_filt)
seurat_filt <- RunPCA(seurat_filt, npcs = 50)
ElbowPlot(seurat_filt, ndims = 50)

dims_use    <- 1:20
seurat_filt <- FindNeighbors(seurat_filt, dims = dims_use)
seurat_filt <- FindClusters(seurat_filt,  resolution = 0.1)
seurat_filt <- RunUMAP(seurat_filt,       dims = dims_use)

DimPlot(seurat_filt, label = TRUE) + ggtitle("After TCR/BCR removal")
saveRDS(seurat_filt, "merged_after_reclustering_no_TCR_BCR.rds")

setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(DESeq2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# Load filtered object
# ============================================================

seurat_filt <- readRDS("merged_after_reclustering_no_TCR_BCR.rds")
cat("Loaded:", nrow(seurat_filt), "genes |", ncol(seurat_filt), "cells\n")

# ============================================================
# Step 1 — FindAllMarkers
# ============================================================

Idents(seurat_filt) <- "RNA_snn_res.0.1"

markers <- FindAllMarkers(
  seurat_filt,
  only.pos        = FALSE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

markers %>%
  filter(p_val_adj < 0.05) %>%
  count(cluster) %>%
  print()

top_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

write.csv(markers, "findallmarkers_res0.1_noTCR_BCR.csv", row.names = FALSE)

# ============================================================
# Step 2 — Setup for pseudobulk DE
# ============================================================

# Recode tissue to avoid underscore in column names
seurat_filt$tissue_clean <- dplyr::recode(seurat_filt$tissue,
  "lymph_node" = "LN",
  "PBMC"       = "PBMC"
)

# Define clusters and maps
clusters_of_interest <- as.character(sort(unique(seurat_filt$RNA_snn_res.0.1)))
sex_map              <- seurat_filt[[]] %>% distinct(patient_id, sex)
study_map            <- seurat_filt[[]] %>% distinct(patient_id, study)

# Verify column name parsing before running loop
sub_test    <- subset(seurat_filt, RNA_snn_res.0.1 == clusters_of_interest[1])
counts_test <- AggregateExpression(sub_test,
  group.by = c("patient_id", "vaccine", "timepoint", "tissue_clean"),
  assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA

cat("Example column names:\n")
print(colnames(counts_test)[1:10])

data.frame(
  colname   = colnames(counts_test)[1:10],
  patient   = gsub("_.*", "", colnames(counts_test)[1:10]),
  vaccine   = gsub(".*_([^_]+)_[^_]+_[^_]+$", "\\1", colnames(counts_test)[1:10]),
  timepoint = gsub(".*_([^_]+)_[^_]+$", "\\1", colnames(counts_test)[1:10]),
  tissue    = gsub(".*_", "", colnames(counts_test)[1:10])
) %>% print()

# ============================================================
# Step 3 — Pseudobulk DE loop
# ============================================================

results_vaccine <- list()
n_clusters      <- length(clusters_of_interest)

for (i in seq_along(clusters_of_interest)) {
  cl <- clusters_of_interest[i]
  cat(sprintf("\n[%d/%d] Processing cluster %s...\n", i, n_clusters, cl))

  sub <- subset(seurat_filt, RNA_snn_res.0.1 == cl)
  cat(sprintf("  > %d cells found\n", ncol(sub)))

  counts_cl <- AggregateExpression(sub,
    group.by = c("patient_id", "vaccine", "timepoint", "tissue_clean"),
    assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA

  meta_cl <- data.frame(
    sample    = colnames(counts_cl),
    patient   = gsub("_.*", "", colnames(counts_cl)),
    vaccine   = gsub(".*_([^_]+)_[^_]+_[^_]+$", "\\1", colnames(counts_cl)),
    timepoint = gsub(".*_([^_]+)_[^_]+$", "\\1", colnames(counts_cl)),
    tissue    = gsub(".*_", "", colnames(counts_cl))
  ) %>%
    mutate(vaccine = factor(vaccine, levels = c("Fluarix", "mRNA-1010", "covidmRNA"))) %>%
    left_join(sex_map,   by = c("patient" = "patient_id")) %>%
    left_join(study_map, by = c("patient" = "patient_id")) %>%
    mutate(
      tissue = factor(tissue), tissue = droplevels(tissue)
    )

  if (any(is.na(meta_cl$tissue))) {
    n_missing <- sum(is.na(meta_cl$tissue))
    cat(sprintf("  > WARNING: %d samples missing tissue — dropping\n", n_missing))
    meta_cl   <- filter(meta_cl, !is.na(tissue)) %>%
      mutate(tissue = droplevels(tissue))
    counts_cl <- counts_cl[, meta_cl$sample]
  }

  # Design: ~ tissue + vaccine, fall back if tissue has <2 levels or is confounded
  has_tissue <- nlevels(meta_cl$tissue) > 1

  tissue_confounded <- if (has_tissue) {
    any(colSums(table(meta_cl$tissue, meta_cl$vaccine) == 0) > 0)
  } else FALSE

  if (!has_tissue) {
    cat("  > NOTE: only one tissue level — using ~ vaccine only\n")
    design_cl <- ~ vaccine
  } else if (tissue_confounded) {
    cat("  > NOTE: tissue confounded with vaccine — using ~ vaccine only\n")
    print(table(meta_cl$tissue, meta_cl$vaccine))
    design_cl <- ~ vaccine
  } else {
    design_cl <- ~ tissue + vaccine
  }

  cat(sprintf("  > Design: %s\n", deparse(design_cl)))

  tryCatch({
    dds_cl <- DESeqDataSetFromMatrix(counts_cl, meta_cl, design = design_cl)
    dds_cl <- DESeq(dds_cl)

    results_vaccine[[cl]] <- list(
      mRNA1010_vs_Fluarix   = results(dds_cl, contrast = c("vaccine", "mRNA-1010", "Fluarix"),    tidy = TRUE),
      covidmRNA_vs_Fluarix  = results(dds_cl, contrast = c("vaccine", "covidmRNA", "Fluarix"),    tidy = TRUE),
      covidmRNA_vs_mRNA1010 = results(dds_cl, contrast = c("vaccine", "covidmRNA", "mRNA-1010"),  tidy = TRUE)
    )

    for (ct in names(results_vaccine[[cl]])) {
      n_sig <- sum(results_vaccine[[cl]][[ct]]$padj < 0.05 &
                   abs(results_vaccine[[cl]][[ct]]$log2FoldChange) > 1,
                   na.rm = TRUE)
      cat(sprintf("     %-30s — %d significant DEGs\n", ct, n_sig))
    }

    cat(sprintf("  [%d/%d] Cluster %s — DONE\n", i, n_clusters, cl))

  }, error = function(e) {
    cat(sprintf("  [%d/%d] Cluster %s — FAILED: %s\n", i, n_clusters, cl, conditionMessage(e)))
  })
}

cat("\n=== All clusters processed ===\n")
cat(sprintf("Successful: %d/%d clusters\n", length(results_vaccine), n_clusters))
cat("Clusters completed:", paste(names(results_vaccine), collapse = ", "), "\n")

# ============================================================
# Step 4 — Flatten results and summarize
# ============================================================

deg_vaccine <- purrr::imap_dfr(results_vaccine, function(cluster_res, cl) {
  purrr::imap_dfr(cluster_res, function(df, contrast_name) {
    df %>% mutate(cluster = cl, contrast = contrast_name)
  })
}) %>%
  rename(gene = row) %>%
  filter(!is.na(padj)) %>%
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1)

deg_vaccine %>%
  filter(significant) %>%
  count(cluster, contrast) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = n, values_fill = 0) %>%
  print()

# ============================================================
# Step 5 — Volcano plots
# ============================================================

plot_volcano <- function(df, cluster_id, contrast_id, top_n = 15) {
  dat <- df %>%
    filter(cluster == cluster_id, contrast == contrast_id) %>%
    mutate(direction = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "up",
      padj < 0.05 & log2FoldChange < -1 ~ "down",
      TRUE ~ "ns"
    ))

  top_genes  <- dat %>% filter(direction != "ns") %>% slice_min(padj, n = top_n) %>% pull(gene)
  dat_labels <- dat %>% filter(gene %in% top_genes)

  ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = dat_labels, aes(label = gene), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c(up = "firebrick", down = "steelblue", ns = "grey70"), drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = paste("Cluster", cluster_id, "|", contrast_id)) +
    theme_bw()
}

pdf("volcano_grids_noTCR_BCR.pdf", width = 18, height = 14)

purrr::walk(unique(deg_vaccine$contrast), function(ct) {
  plots <- purrr::map(clusters_of_interest, function(cl) {
    dat <- deg_vaccine %>% filter(cluster == cl, contrast == ct)

    if (nrow(dat) == 0) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5,
          label = paste("Cluster", cl, "\nno data"),
          hjust = 0.5, vjust = 0.5, size = 4, color = "grey50") +
        theme_void())
    }

    plot_volcano(deg_vaccine, cl, ct) +
      labs(title = paste("Cluster", cl)) +
      theme(legend.position = "none") +
      scale_color_manual(
        values = c(up = "firebrick", down = "steelblue", ns = "grey70"),
        drop   = FALSE)
  })

  p_grid <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = ct,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  print(p_grid)
})

dev.off()
cat("Saved to volcano_grids_noTCR_BCR.pdf\n")
