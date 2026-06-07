setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

# ============================================================
# Pseudobulk DE + KEGG + MSigDB C2 enrichment
# Overall DE and per-cluster DE
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

# ─────────────────────────────────────────────
# Load object and add cell-type annotation
# ─────────────────────────────────────────────

seurat_obj <- readRDS("merged_after_reclustering_no_TCR_BCR.rds")

cluster_to_celltype <- c(
  "0"  = "Naive/central memory CD4 T",
  "1"  = "Activated/memory B cells",
  "2"  = "Naive B",
  "3"  = "Plasma cells",
  "4"  = "Cytotoxic CD8 T",
  "5"  = "Germinal center B",
  "6"  = "Cytotoxic lymphocytes", # mixed NK and cytotoxic T at res 0.1
  "7"  = "T follicular helper (Tfh)",
  "8"  = "Classical monocytes",
  "9"  = "Activated/memory T",
  "10" = "Plasmablasts",
  "11" = "Plasmacytoid dendritic cells (pDC)"
)

cell_type <- cluster_to_celltype[as.character(seurat_obj$RNA_snn_res.0.1)]
names(cell_type) <- colnames(seurat_obj)

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = cell_type,
  col.name = "cell_type"
)

print(table(seurat_obj$cell_type, useNA = "ifany"))

all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.25
)

top3_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3) %>%
  ungroup()

marker_genes <- unique(top3_markers$gene)

marker_genes <- top3_markers %>%
  filter(!grepl("^LINC", gene)) %>%
  arrange(as.numeric(cluster), desc(avg_log2FC)) %>%
  pull(gene)

DotPlot(
  seurat_obj,
  features = marker_genes,
  group.by = "RNA_snn_res.0.1"
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

cluster_labels <- seurat_obj@meta.data %>%
  dplyr::select(RNA_snn_res.0.1, cell_type) %>%
  distinct() %>%
  arrange(as.numeric(RNA_snn_res.0.1))

cluster_labels

cluster_names <- setNames(
  cluster_labels$cell_type,
  cluster_labels$RNA_snn_res.0.1
)

seurat_obj$cluster_label <- factor(
  seurat_obj$RNA_snn_res.0.1,
  levels = names(cluster_names),
  labels = cluster_names
)

DotPlot(
  seurat_obj,
  features = marker_genes,
  group.by = "cluster_label"
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )
# ─────────────────────────────────────────────
# Basic UMAP checks
# ─────────────────────────────────────────────

DimPlot(seurat_obj, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "vaccine")

# ─────────────────────────────────────────────
# Composition plots
# ─────────────────────────────────────────────

meta <- seurat_obj[[]]

df <- seurat_obj@meta.data %>%
  select(cluster = RNA_snn_res.0.1, study, vaccine, sex, tissue)

df_long <- bind_rows(
  df %>% count(cluster, group = study)   %>% mutate(variable = "study"),
  df %>% count(cluster, group = vaccine) %>% mutate(variable = "vaccine"),
  df %>% count(cluster, group = sex)     %>% mutate(variable = "sex"),
  df %>% count(cluster, group = tissue)  %>% mutate(variable = "tissue")
) %>%
  group_by(variable, cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_long, aes(x = cluster, y = prop, fill = group)) +
  geom_col(position = "fill") +
  facet_wrap(~ variable) +
  ylab("Proportion within cluster") +
  xlab("Cluster, resolution 0.1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_pseudobulk_composition <- function(meta,
                                        sample_col  = "patient_id",
                                        cluster_col = "RNA_snn_res.0.1",
                                        group_col   = NULL,
                                        min_cells   = 10) {
  
  grouping_vars <- c(sample_col, cluster_col)
  if (!is.null(group_col)) grouping_vars <- c(grouping_vars, group_col)
  if ("tissue" %in% colnames(meta) && !"tissue" %in% grouping_vars) {
    grouping_vars <- c(grouping_vars, "tissue")
  }
  
  counts <- meta %>%
    filter(!is.na(.data[[sample_col]])) %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  props <- counts %>%
    group_by(across(all_of(setdiff(grouping_vars, cluster_col)))) %>%
    mutate(total = sum(n_cells), proportion = n_cells / total) %>%
    filter(total >= min_cells) %>%
    ungroup()
  
  p <- ggplot(props, aes(x = .data[[cluster_col]], y = proportion)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(aes(color = .data[[sample_col]]), width = 0.15, size = 2, alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    labs(
      x = "Cluster",
      y = "Proportion of cells (per donor × tissue)",
      color = sample_col,
      title = paste("Pseudobulk cluster composition",
                    if (!is.null(group_col)) paste("by", group_col) else "")
    )
  
  if (!is.null(group_col) && "tissue" %in% colnames(meta)) {
    p <- p + facet_grid(rows = vars(.data[[group_col]]), cols = vars(tissue))
  } else if ("tissue" %in% colnames(meta)) {
    p <- p + facet_wrap(~ tissue)
  } else if (!is.null(group_col)) {
    p <- p + facet_wrap(~ .data[[group_col]])
  }
  
  p
}

plot_pseudobulk_composition(meta, group_col = "vaccine")
plot_pseudobulk_composition(meta, group_col = "timepoint")
plot_pseudobulk_composition(meta, group_col = "study")
plot_pseudobulk_composition(meta, group_col = "sex")


study_table <- meta %>%
  distinct(patient_id, study, tissue, timepoint, vaccine) %>%
  count(study, vaccine, tissue, timepoint, name = "n_donors") %>%
  arrange(study, vaccine,  tissue, timepoint)

View(study_table)


# ─────────────────────────────────────────────
# Helper: pseudobulk aggregation
# ─────────────────────────────────────────────

make_pseudobulk <- function(seurat_obj,
                            group_vars,
                            assay = "RNA") {
  
  counts <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  
  meta <- seurat_obj[[]] %>%
    filter(!is.na(patient_id), !is.na(vaccine))
  
  valid_cells <- rownames(meta)[complete.cases(meta[, group_vars, drop = FALSE])]
  meta <- meta[valid_cells, , drop = FALSE]
  counts <- counts[, valid_cells, drop = FALSE]
  
  meta$pb_group <- apply(meta[, group_vars, drop = FALSE], 1, paste, collapse = "__")
  
  groups <- unique(meta$pb_group)
  
  pb_counts <- sapply(groups, function(g) {
    cells <- rownames(meta)[meta$pb_group == g]
    if (length(cells) == 1) {
      as.numeric(counts[, cells])
    } else {
      Matrix::rowSums(counts[, cells, drop = FALSE])
    }
  })
  
  rownames(pb_counts) <- rownames(counts)
  colnames(pb_counts) <- groups
  
  pb_meta <- meta %>%
    rownames_to_column("cell") %>%
    group_by(pb_group) %>%
    slice(1) %>%
    ungroup() %>%
    column_to_rownames("pb_group") %>%
    select(all_of(group_vars))
  
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  
  list(counts = pb_counts, meta = pb_meta)
}

# ─────────────────────────────────────────────
# Helper: clean canonical contrasts
# Keeps:
#   covidmRNA vs Fluarix
#   mRNA-1010 vs Fluarix
#   mRNA-1010 vs covidmRNA
# Drops:
#   Fluarix vs covidmRNA
#   Fluarix vs mRNA-1010
# ─────────────────────────────────────────────

clean_contrasts <- function(de_tbl) {
  de_tbl %>%
    distinct(gene, tissue, across(any_of("cluster")), contrast, .keep_all = TRUE) %>%
    filter(!contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010"))
}

# ─────────────────────────────────────────────
# Overall DE
# ─────────────────────────────────────────────

run_overall_de <- function(seurat_obj,
                           contrast_var = "vaccine",
                           ref_level    = "Fluarix",
                           min_counts   = 10) {
  
  message("Creating overall pseudobulk object...")
  
  pb <- make_pseudobulk(
    seurat_obj,
    group_vars = c("patient_id", "tissue", "vaccine")
  )
  
  results_list <- list()
  
  for (tis in sort(unique(pb$meta$tissue))) {
    
    message("\n========================================")
    message("Overall DE; tissue: ", tis)
    
    keep <- pb$meta$tissue == tis
    
    counts <- pb$counts[, keep, drop = FALSE]
    s_meta <- pb$meta[keep, , drop = FALSE]
    
    keep_samples <- colSums(counts) >= min_counts
    counts <- counts[, keep_samples, drop = FALSE]
    s_meta <- s_meta[keep_samples, , drop = FALSE]
    
    if (ncol(counts) < 4) {
      message("Skipping ", tis, ": <4 samples")
      next
    }
    
    levels_present <- unique(s_meta[[contrast_var]])
    
    if (!ref_level %in% levels_present) {
      message("Skipping ", tis, ": reference not present")
      next
    }
    
    if (length(levels_present) < 2) {
      message("Skipping ", tis, ": only one vaccine level")
      next
    }
    
    s_meta[[contrast_var]] <- relevel(factor(s_meta[[contrast_var]]), ref = ref_level)
    
    keep_genes <- rowSums(counts >= 1) >= 2
    counts <- counts[keep_genes, , drop = FALSE]
    
    if (nrow(counts) < 10) {
      message("Skipping ", tis, ": too few genes")
      next
    }
    
    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = round(counts),
        colData   = s_meta,
        design    = as.formula(paste("~", contrast_var))
      )
      
      dds <- DESeq(dds, quiet = TRUE)
      
      for (lvl in setdiff(levels(s_meta[[contrast_var]]), ref_level)) {
        
        res <- results(dds, contrast = c(contrast_var, lvl, ref_level)) %>%
          as.data.frame() %>%
          rownames_to_column("gene") %>%
          mutate(
            tissue = tis,
            contrast = paste(lvl, "vs", ref_level)
          )
        
        results_list[[paste(tis, lvl, ref_level, sep = "__")]] <- res
      }
    }, error = function(e) {
      message("Error in tissue=", tis, ": ", e$message)
    })
  }
  
  bind_rows(results_list)
}

de_overall_fluarix_ref <- run_overall_de(seurat_obj, ref_level = "Fluarix")
de_overall_covid_ref   <- run_overall_de(seurat_obj, ref_level = "covidmRNA")

de_overall_all <- bind_rows(de_overall_fluarix_ref, de_overall_covid_ref) %>%
  clean_contrasts()

#saveRDS(de_overall_all, "pseudobulk_DE_overall_all_contrasts.rds")
write.csv(de_overall_all, "pseudobulk_DE_overall_all_contrasts.csv", row.names = FALSE)

overall_summary_tbl <- de_overall_all %>%
  filter(!is.na(padj)) %>%
  group_by(tissue, contrast) %>%
  summarise(
    n_sig_up   = sum(padj < 0.05 & log2FoldChange >  0.5),
    n_sig_down = sum(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  ) %>%
  arrange(tissue, contrast)

write.csv(overall_summary_tbl, "de_overall_summary.csv", row.names = FALSE)

# ─────────────────────────────────────────────
# Overall DE volcano plot
# ─────────────────────────────────────────────

de_plot <- de_overall_all %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  0.5 ~ "up",
    padj < 0.05 & log2FoldChange < -0.5 ~ "down",
    TRUE ~ "ns"
  ))

top_genes <- de_plot %>%
  filter(sig != "ns") %>%
  group_by(tissue, contrast) %>%
  slice_min(padj, n = 10) %>%
  ungroup()

p_overall_volcano <- ggplot(de_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 2.2,
    max.overlaps = 15,
    segment.size = 0.3,
    color = "black"
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(up = "#e63946", down = "#457b9d", ns = "grey75")) +
  scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
  facet_grid(rows = vars(tissue), cols = vars(contrast)) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = NULL,
    title = "Pseudobulk DE — vaccine comparisons by tissue"
  )

ggsave("DE_overall_volcano_by_tissue.pdf", p_overall_volcano, width = 14, height = 8)

# ─────────────────────────────────────────────
# KEGG enrichment helper and overall KEGG
# ─────────────────────────────────────────────

run_kegg <- function(genes, universe = NULL, pval = 0.05) {
  
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  
  if (length(genes) < 5) return(NULL)
  
  entrez <- suppressMessages(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  
  if (is.null(entrez) || nrow(entrez) < 5) return(NULL)
  
  uni_entrez <- NULL
  if (!is.null(universe)) {
    uni_entrez <- suppressMessages(
      bitr(unique(universe), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )$ENTREZID
  }
  
  enr <- enrichKEGG(
    gene = unique(entrez$ENTREZID),
    universe = unique(uni_entrez),
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = pval,
    qvalueCutoff = 0.2
  )
  
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  enr
}

run_enrichment_by_group <- function(de_tbl,
                                    run_fun,
                                    output_prefix,
                                    group_cols,
                                    universe = NULL,
                                    padj_cutoff = 0.05,
                                    lfc_cutoff = 0.5,
                                    top_n = 5) {
  
  de_for_enrich <- de_tbl %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      padj < padj_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    ) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down"))
  
  if (nrow(de_for_enrich) == 0) {
    warning("No significant genes for enrichment: ", output_prefix)
    return(NULL)
  }
  
  results <- de_for_enrich %>%
    group_by(across(all_of(c(group_cols, "direction")))) %>%
    group_map(~ {
      res <- run_fun(.x$gene, universe = universe)
      if (is.null(res)) return(NULL)
      
      as.data.frame(res) %>%
        mutate(
          !!!as.list(.y),
          gene_input_n = length(unique(.x$gene))
        )
    }, .keep = TRUE) %>%
    bind_rows()
  
  if (is.null(results) || nrow(results) == 0) {
    warning("No enriched terms returned for: ", output_prefix)
    return(NULL)
  }
  
  write.csv(results, paste0(output_prefix, "_enrichment.csv"), row.names = FALSE)
  saveRDS(results, paste0(output_prefix, "_enrichment.rds"))
  
  top_paths <- results %>%
    group_by(across(all_of(c(group_cols, "direction")))) %>%
    slice_min(p.adjust, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  p <- ggplot(
    top_paths,
    aes(
      x = direction,
      y = reorder(Description, -p.adjust),
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point() +
    scale_color_gradient(low = "#e63946", high = "#adb5bd", name = "adj. p-value") +
    scale_size_continuous(range = c(2, 8), name = "Gene count") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 7),
      strip.background = element_rect(fill = "grey92"),
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "right"
    ) +
    labs(x = "Direction", y = NULL, title = paste0(output_prefix, " enrichment"))
  
  if (identical(group_cols, c("tissue", "contrast"))) {
    p <- p + facet_grid(
      cols = vars(tissue),
      rows = vars(contrast),
      scales = "free_y",
      space = "free_y"
    )
  } else if (identical(group_cols, c("tissue", "cluster", "contrast"))) {
    p <- p + facet_grid(
      rows = vars(contrast),
      cols = vars(tissue),
      scales = "free_y",
      space = "free_y"
    )
  }
  
  ggsave(paste0(output_prefix, "_dotplot.pdf"), p, width = 14, height = 10)
  
  list(results = results, plot = p)
}

kegg_overall <- run_enrichment_by_group(
  de_tbl = de_overall_all,
  run_fun = run_kegg,
  output_prefix = "kegg_overall",
  group_cols = c("tissue", "contrast"),
  universe = unique(de_overall_all$gene)
)

# ─────────────────────────────────────────────
# Per-cluster DE
# ─────────────────────────────────────────────

run_percluster_de <- function(seurat_obj,
                              cluster_col  = "RNA_snn_res.0.1",
                              contrast_var = "vaccine",
                              ref_level    = "Fluarix",
                              min_counts   = 10) {
  
  message("Creating per-cluster pseudobulk object...")
  
  seurat_obj$pb_cluster <- seurat_obj[[cluster_col]][, 1]
  
  pb <- make_pseudobulk(
    seurat_obj,
    group_vars = c("patient_id", "tissue", "vaccine", "pb_cluster")
  )
  
  results_list <- list()
  
  for (tis in sort(unique(pb$meta$tissue))) {
    
    message("\n========================================")
    message("Per-cluster DE; tissue: ", tis)
    
    for (cl in sort(unique(pb$meta$pb_cluster))) {
      
      message("----------------------------------------")
      message("Cluster: ", cl)
      
      keep <- pb$meta$tissue == tis & pb$meta$pb_cluster == cl
      
      counts <- pb$counts[, keep, drop = FALSE]
      s_meta <- pb$meta[keep, , drop = FALSE]
      
      if (ncol(counts) < 4) {
        message("Skipping: <4 pseudobulk samples")
        next
      }
      
      keep_samples <- colSums(counts) >= min_counts
      counts <- counts[, keep_samples, drop = FALSE]
      s_meta <- s_meta[keep_samples, , drop = FALSE]
      
      if (ncol(counts) < 4) {
        message("Skipping after count filter: <4 pseudobulk samples")
        next
      }
      
      levels_present <- unique(s_meta[[contrast_var]])
      
      if (!ref_level %in% levels_present) {
        message("Skipping: reference level not present")
        next
      }
      
      if (length(levels_present) < 2) {
        message("Skipping: only one vaccine level")
        next
      }
      
      s_meta[[contrast_var]] <- relevel(factor(s_meta[[contrast_var]]), ref = ref_level)
      
      keep_genes <- rowSums(counts >= 1) >= 2
      counts <- counts[keep_genes, , drop = FALSE]
      
      if (nrow(counts) < 10) {
        message("Skipping: too few genes")
        next
      }
      
      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = round(counts),
          colData   = s_meta,
          design    = as.formula(paste("~", contrast_var))
        )
        
        dds <- DESeq(dds, quiet = TRUE)
        
        for (lvl in setdiff(levels(s_meta[[contrast_var]]), ref_level)) {
          
          res <- results(dds, contrast = c(contrast_var, lvl, ref_level)) %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(
              tissue = tis,
              cluster = cl,
              cell_type = unname(cluster_to_celltype[as.character(cl)]),
              contrast = paste(lvl, "vs", ref_level)
            )
          
          results_list[[paste(tis, cl, lvl, ref_level, sep = "__")]] <- res
        }
      }, error = function(e) {
        message("Error in tissue=", tis, " cluster=", cl, ": ", e$message)
      })
    }
  }
  
  bind_rows(results_list)
}

de_percluster_fluarix_ref <- run_percluster_de(seurat_obj, ref_level = "Fluarix")
de_percluster_covid_ref   <- run_percluster_de(seurat_obj, ref_level = "covidmRNA")

de_percluster_all <- bind_rows(de_percluster_fluarix_ref, de_percluster_covid_ref) %>%
  clean_contrasts()

#saveRDS(de_percluster_all, "pseudobulk_DE_percluster_all_contrasts.rds")
write.csv(de_percluster_all, "pseudobulk_DE_percluster_all_contrasts.csv", row.names = FALSE)

summary_tbl <- de_percluster_all %>%
  filter(!is.na(padj)) %>%
  group_by(tissue, cluster, cell_type, contrast) %>%
  summarise(
    n_sig_up   = sum(padj < 0.05 & log2FoldChange >  0.5),
    n_sig_down = sum(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  ) %>%
  arrange(tissue, cluster, contrast)

write.csv(summary_tbl, "de_percluster_summary.csv", row.names = FALSE)

# ─────────────────────────────────────────────
# Per-cluster volcano PDF
# ─────────────────────────────────────────────

de_sig <- de_percluster_all %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  0.5 ~ "up",
    padj < 0.05 & log2FoldChange < -0.5 ~ "down",
    TRUE ~ "ns"
  ))

top_genes_pc <- de_sig %>%
  filter(sig != "ns") %>%
  group_by(tissue, contrast, cluster) %>%
  slice_min(padj, n = 10) %>%
  ungroup()

cairo_pdf("DE_percluster_volcanos_by_contrast.pdf", width = 20, height = 16, onefile = TRUE)

for (ct in unique(de_sig$contrast)) {
  for (tis in sort(unique(de_sig$tissue))) {
    
    dat <- filter(de_sig, contrast == ct, tissue == tis)
    dat_labs <- filter(top_genes_pc, contrast == ct, tissue == tis)
    
    if (nrow(dat) == 0) next
    
    counts_lab <- dat %>%
      filter(sig != "ns") %>%
      count(cluster, sig) %>%
      pivot_wider(names_from = sig, values_from = n, values_fill = 0) %>%
      mutate(
        up = ifelse("up" %in% colnames(.), up, 0),
        down = ifelse("down" %in% colnames(.), down, 0),
        label = paste0("\u2191", up, "  \u2193", down)
      )
    
    p <- ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text_repel(
        data = dat_labs,
        aes(label = gene),
        size = 2.2,
        max.overlaps = 12,
        segment.size = 0.3,
        color = "black"
      ) +
      geom_vline(xintercept = c(-0.5, 0.5),
                 linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_text(
        data = counts_lab,
        aes(label = label, x = Inf, y = Inf),
        hjust = 1.1,
        vjust = 1.5,
        size = 3,
        color = "grey30",
        inherit.aes = FALSE
      ) +
      scale_color_manual(values = c(up = "#e63946", down = "#457b9d", ns = "grey75")) +
      scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
      facet_wrap(~ cluster, ncol = 4, labeller = label_both) +
      theme_bw(base_size = 11) +
      theme(
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40")
      ) +
      labs(
        x = "log2 Fold Change",
        y = "-log10(adjusted p-value)",
        color = NULL,
        title = paste0("Contrast: ", ct),
        subtitle = paste0("Tissue: ", tis, " — one panel per cluster")
      )
    
    print(p)
  }
}

dev.off()
message("Saved: C2_DE_percluster_volcanos_by_contrast.pdf")

# ============================================================
# MSigDB C2 enrichment
# Runs for:
#   1. Overall DEGs
#   2. Per-cluster DEGs
# ============================================================

# The current msigdbr API uses collection/subcollection.
# This fallback keeps the script compatible with older msigdbr versions.
get_msigdb_c2 <- function() {
  args <- names(formals(msigdbr::msigdbr))
  
  if ("collection" %in% args) {
    msigdbr(species = "Homo sapiens", collection = "C2")
  } else {
    msigdbr(species = "Homo sapiens", category = "C2")
  }
}

msigdb_c2 <- get_msigdb_c2() %>%
  select(gs_name, gene_symbol) %>%
  distinct()

run_msigdb_c2 <- function(genes, universe = NULL, pval = 0.05) {
  
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  
  if (length(genes) < 5) return(NULL)
  
  if (!is.null(universe)) {
    universe <- unique(universe)
    universe <- universe[!is.na(universe)]
  }
  
  enr <- enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = msigdb_c2,
    pAdjustMethod = "BH",
    pvalueCutoff = pval,
    qvalueCutoff = 0.2
  )
  
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  enr
}

# ─────────────────────────────────────────────
# MSigDB C2: overall DEGs
# ─────────────────────────────────────────────

msigdb_c2_overall_obj <- run_enrichment_by_group(
  de_tbl = de_overall_all,
  run_fun = run_msigdb_c2,
  output_prefix = "msigdb_c2_overall",
  group_cols = c("tissue", "contrast"),
  universe = unique(de_overall_all$gene),
  padj_cutoff = 0.05,
  lfc_cutoff = 0.5,
  top_n = 5
)

msigdb_c2_overall <- if (!is.null(msigdb_c2_overall_obj)) msigdb_c2_overall_obj$results else NULL

# ─────────────────────────────────────────────
# MSigDB C2: per-cluster DEGs
# This is the new part added for cluster-level enrichment.
# It runs separately by tissue × cluster × contrast × direction.
# ─────────────────────────────────────────────

msigdb_c2_percluster_obj <- run_enrichment_by_group(
  de_tbl = de_percluster_all,
  run_fun = run_msigdb_c2,
  output_prefix = "msigdb_c2_percluster",
  group_cols = c("tissue", "cluster", "contrast"),
  universe = unique(de_percluster_all$gene),
  padj_cutoff = 0.05,
  lfc_cutoff = 0.5,
  top_n = 5
)

msigdb_c2_percluster <- if (!is.null(msigdb_c2_percluster_obj)) msigdb_c2_percluster_obj$results else NULL

# Optional cleaner per-cluster C2 plot: one PDF page per tissue/contrast.
# This is often more readable than trying to facet all clusters in one panel.

if (!is.null(msigdb_c2_percluster) && nrow(msigdb_c2_percluster) > 0) {
  
  cairo_pdf("msigdb_c2_percluster_dotplots_by_tissue_contrast.pdf",
            width = 18, height = 12, onefile = TRUE)
  
  for (tis in sort(unique(msigdb_c2_percluster$tissue))) {
    for (ct in unique(msigdb_c2_percluster$contrast)) {
      
      dat <- msigdb_c2_percluster %>%
        filter(tissue == tis, contrast == ct) %>%
        group_by(cluster, direction) %>%
        slice_min(p.adjust, n = 5, with_ties = FALSE) %>%
        ungroup()
      
      if (nrow(dat) == 0) next
      
      p <- ggplot(
        dat,
        aes(
          x = direction,
          y = reorder(Description, -p.adjust),
          size = Count,
          color = p.adjust
        )
      ) +
        geom_point() +
        scale_color_gradient(low = "#e63946", high = "#adb5bd", name = "adj. p-value") +
        scale_size_continuous(range = c(2, 8), name = "Gene count") +
        facet_wrap(~ cluster, scales = "free_y", ncol = 3, labeller = label_both) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.y = element_text(size = 6),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold", size = 8),
          legend.position = "right",
          panel.grid.minor = element_blank()
        ) +
        labs(
          x = "Direction",
          y = NULL,
          title = paste0("MSigDB C2 enrichment — ", ct),
          subtitle = paste0("Tissue: ", tis, " — one panel per cluster")
        )
      
      print(p)
    }
  }
  
  dev.off()
  message("Saved: msigdb_c2_percluster_dotplots_by_tissue_contrast.pdf")
}

# ─────────────────────────────────────────────
# Save final workspace objects
# ─────────────────────────────────────────────

saveRDS(de_overall_all, "pseudobulk_DE_overall.rds")
saveRDS(de_percluster_all, "pseudobulk_DE_percluster.rds")
saveRDS(msigdb_c2_overall, "msigdb_c2_overall_enrichment.rds")
saveRDS(msigdb_c2_percluster, "msigdb_c2_percluster_enrichment.rds")

message("Done.")
