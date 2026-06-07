setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(DESeq2)
library(volcano3D)
library(Seurat)
library(plotly)
library(DESeq2)

####adding metadata

seurat_obj <- readRDS("merged_after_reclustering_no_TCR_BCR.rds")

cluster_to_celltype <- c(
  "0"  = "Naive/central memory CD4 T",
  "1"  = "Activated/memory B cells",
  "2"  = "Naive B",
  "3"  = "Plasma cells",
  "4"  = "Cytotoxic CD8 T",
  "5"  = "Germinal center B",
  "6"  = "Cytotoxic Lymphocites", #mix of nk and cytotoxic t cells not distinguished in res 0.1
  "7"  = "T follicular helper (Tfh)",
  "8"  = "Classical monocytes",
  "9"  = "Activated/memory T",
  "10" = "Plasmablasts",
  "11" = "Plasmacytoid dendritic cells (pDC)"
)

cell_type <- cluster_to_celltype[
  as.character(seurat_obj$RNA_snn_res.0.1)
]

names(cell_type) <- colnames(seurat_obj)

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = cell_type,
  col.name = "cell_type"
)

table(seurat_obj$cell_type, useNA = "ifany")

FindMarkers(
  seurat_obj,
  ident.1 = "1",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
) %>%
  arrange(desc(avg_log2FC)) %>%
  head(50)

meta <- seurat_obj[[]]

DimPlot(seurat_obj, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "vaccine")

library(dplyr)
library(tidyr)
library(ggplot2)

df <- seurat_obj@meta.data %>%
  select(cluster = RNA_snn_res.0.1, study, vaccine, sex, tissue)

df_long <- bind_rows(
  df %>% count(cluster, group = study) %>% mutate(variable = "study"),
  df %>% count(cluster, group = vaccine) %>% mutate(variable = "vaccine"),
  df %>% count(cluster, group = sex) %>% mutate(variable = "sex"),
  df %>% count(cluster, group =  tissue) %>% mutate(variable = "tissue")
  ) %>%
  group_by(variable, cluster) %>%
  mutate(prop = n / sum(n))

ggplot(df_long, aes(x = cluster, y = prop, fill = group)) +
  geom_col(position = "fill") +
  facet_wrap(~variable) +
  ylab("Proportion within cluster") +
  xlab("Cluster, resolution 0.1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#############pseudobulk Compositional plots
library(dplyr)
library(ggplot2)
library(tidyr)

plot_pseudobulk_composition <- function(meta,
                                        sample_col   = "patient_id",
                                        cluster_col  = "RNA_snn_res.0.1",
                                        group_col    = NULL,
                                        min_cells    = 10) {
  
  # 1. Count cells per sample × tissue × cluster
  #    Always include tissue in the grouping if it exists
  grouping_vars <- c(sample_col, cluster_col)
  if (!is.null(group_col)) grouping_vars <- c(grouping_vars, group_col)
  if ("tissue" %in% colnames(meta) && !"tissue" %in% grouping_vars) {
    grouping_vars <- c(grouping_vars, "tissue")
  }
  
  counts <- meta %>%
    filter(!is.na(.data[[sample_col]])) %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(n_cells = n(), .groups = "drop")
  
  # 2. Proportion within each sample × tissue
  props <- counts %>%
    group_by(across(all_of(setdiff(grouping_vars, cluster_col)))) %>%
    mutate(
      total      = sum(n_cells),
      proportion = n_cells / total
    ) %>%
    filter(total >= min_cells) %>%
    ungroup()
  
  # 3. Plot
  p <- ggplot(props, aes(x = .data[[cluster_col]], y = proportion)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(
      aes(color = .data[[sample_col]]),
      width = 0.15, size = 2, alpha = 0.8
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    labs(
      x     = "Cluster",
      y     = "Proportion of cells (per donor × tissue)",
      color = sample_col,
      title = paste("Pseudobulk cluster composition",
                    if (!is.null(group_col)) paste("by", group_col) else "")
    )
  
  # 4. Facet — tissue always, group_col on top if provided
  if (!is.null(group_col) && "tissue" %in% colnames(meta)) {
    p <- p + facet_grid(rows = vars(.data[[group_col]]),
                        cols = vars(tissue))
  } else if ("tissue" %in% colnames(meta)) {
    p <- p + facet_wrap(~ tissue)
  } else if (!is.null(group_col)) {
    p <- p + facet_wrap(~ .data[[group_col]])
  }
  
  return(p)
}

# Usage examples
plot_pseudobulk_composition(meta, group_col = "vaccine")
plot_pseudobulk_composition(meta, group_col = "timepoint")
plot_pseudobulk_composition(meta, group_col = "study")
plot_pseudobulk_composition(meta, group_col = "sex")

#############################
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(purrr)

# ─────────────────────────────────────────────
# STEP 1: Build pseudobulk count matrix
# ─────────────────────────────────────────────
# Aggregates raw counts per donor × tissue (× cluster for per-cluster analysis)

make_pseudobulk <- function(seurat_obj,
                            group_vars,        # e.g. c("patient_id", "tissue") or add "cluster"
                            assay = "RNA") {
  
  counts <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  meta   <- seurat_obj[[]] %>%
    filter(!is.na(patient_id), !is.na(vaccine))
  
  # Build a grouping key per cell
  meta$pb_group <- apply(meta[, group_vars, drop = FALSE], 1, paste, collapse = "__")
  
  # Keep only cells with all grouping vars populated
  valid_cells <- rownames(meta)[complete.cases(meta[, group_vars])]
  counts <- counts[, valid_cells]
  meta   <- meta[valid_cells, ]
  
  # Aggregate counts by group
  groups     <- unique(meta$pb_group)
  pb_counts  <- sapply(groups, function(g) {
    cells <- rownames(meta)[meta$pb_group == g]
    if (length(cells) == 1) counts[, cells]
    else Matrix::rowSums(counts[, cells])
  })
  
  # Build sample-level metadata
  pb_meta <- meta %>%
    tibble::rownames_to_column("cell") %>%
    group_by(pb_group) %>%
    slice(1) %>%
    ungroup() %>%
    column_to_rownames("pb_group") %>%
    select(all_of(group_vars))
  
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  
  list(counts = pb_counts, meta = pb_meta)
}

# ─────────────────────────────────────────────
# STEP 2: Overall DE — vaccine comparison
#          stratified by tissue
# ─────────────────────────────────────────────
run_overall_de <- function(seurat_obj,
                           contrast_var  = "vaccine",
                           ref_level     = "Fluarix",
                           min_cells     = 10,
                           min_counts    = 10) {
  
  message("Creating pseudobulk object...")
  pb <- make_pseudobulk(
    seurat_obj,
    group_vars = c("patient_id", "tissue", "vaccine")
  )
  
  message("Pseudobulk created. ",
          ncol(pb$counts), " samples, ",
          nrow(pb$counts), " genes.")
  
  results_list <- list()
  
  for (tis in unique(pb$meta$tissue)) {
    
    message("\n========================================")
    message("Processing tissue: ", tis)
    
    keep <- pb$meta$tissue == tis
    
    counts <- pb$counts[, keep]
    s_meta <- pb$meta[keep, , drop = FALSE]
    
    message("Initial samples: ", ncol(counts))
    
    # Drop low-count pseudobulk samples
    keep2 <- colSums(counts) >= min_counts
    
    counts <- counts[, keep2]
    s_meta <- s_meta[keep2, , drop = FALSE]
    
    message("Samples after count filter (>= ",
            min_counts, "): ", ncol(counts))
    
    # Need ≥ 2 levels of contrast_var
    n_levels <- length(unique(s_meta[[contrast_var]]))
    
    if (n_levels < 2) {
      message("Skipping ", tis,
              ": only ", n_levels,
              " level(s) of ", contrast_var)
      next
    }
    
    message("Levels present: ",
            paste(unique(s_meta[[contrast_var]]), collapse = ", "))
    
    s_meta[[contrast_var]] <- relevel(
      factor(s_meta[[contrast_var]]),
      ref = ref_level
    )
    
    # Filter low-expressed genes
    keep_genes <- rowSums(counts >= 1) >= 2
    
    message("Genes before filtering: ", nrow(counts))
    
    counts <- counts[keep_genes, ]
    
    message("Genes after filtering: ", nrow(counts))
    
    message("Building DESeq2 dataset...")
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData   = s_meta,
      design    = as.formula(paste("~", contrast_var))
    )
    
    message("Running DESeq2...")
    
    dds <- DESeq(dds, quiet = TRUE)
    
    message("DESeq2 finished for tissue: ", tis)
    
    # Extract all pairwise results vs reference
    for (lvl in setdiff(levels(s_meta[[contrast_var]]), ref_level)) {
      
      message("  Extracting contrast: ",
              lvl, " vs ", ref_level)
      
      res <- results(
        dds,
        contrast = c(contrast_var, lvl, ref_level)
      ) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        mutate(
          tissue = tis,
          contrast = paste(lvl, "vs", ref_level)
        )
      
      message("    Retrieved ", nrow(res), " genes")
      
      results_list[[paste(tis, lvl, sep = "__")]] <- res
    }
    
    message("Completed tissue: ", tis)
  }
  
  message("\nCombining results...")
  
  final_results <- bind_rows(results_list)
  
  message("Done. Total rows: ", nrow(final_results))
  
  return(final_results)
}

# Run it
de_overall <- run_overall_de(seurat_obj, ref_level = "Fluarix")

# Run with covidmRNA as reference to get mRNA-1010 vs covidmRNA
de_overall_covid_ref <- run_overall_de(seurat_obj, ref_level = "covidmRNA")

# Bind with original results
de_overall_all <- bind_rows(de_overall, de_overall_covid_ref)

# Now the summary has all three contrasts
de_overall_all %>%
  filter(!is.na(padj)) %>%
  group_by(tissue, contrast) %>%
  summarise(
    n_sig_up   = sum(padj < 0.05 & log2FoldChange > 0.5),
    n_sig_down = sum(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  ) %>%
  mutate(
    contrast_canonical = map_chr(
      str_split(contrast, " vs "),
      ~ paste(sort(.x), collapse = " vs ")
    )
  ) %>%
  distinct(tissue, contrast_canonical, .keep_all = TRUE) %>%
  select(-contrast_canonical) %>%
  arrange(tissue, contrast)

library(ggplot2)
library(dplyr)
library(ggrepel)

# Remove duplicate mirrored contrasts — keep canonical direction only
de_plot <- de_overall_all %>%
  filter(!contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010")) %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# Label top genes per panel
top_genes <- de_plot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  group_by(tissue, contrast) %>%
  slice_min(padj, n = 10) %>%
  ungroup()

# Significance category
de_plot <- de_plot %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  0.5 ~ "up",
    padj < 0.05 & log2FoldChange < -0.5 ~ "down",
    TRUE                                 ~ "ns"
  ))

ggplot(de_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 2.2, max.overlaps = 15,
    segment.size = 0.3, color = "black"
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_color_manual(
    values = c(up = "#e63946", down = "#457b9d", ns = "grey75"),
    labels = c(up = "Up", down = "Down", ns = "NS")
  ) +
  scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
  facet_grid(rows = vars(tissue), cols = vars(contrast)) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x     = "log2 Fold Change",
    y     = "-log10(adjusted p-value)",
    color = NULL,
    title = "Pseudobulk DE — vaccine comparisons by tissue"
  )

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(purrr)

# ─────────────────────────────────────────────
# STEP 1: Helper — run KEGG for one gene list
# ─────────────────────────────────────────────
run_kegg <- function(genes, universe = NULL, pval = 0.05) {
  
  # Convert gene symbols → Entrez IDs
  entrez <- bitr(genes,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
  
  if (nrow(entrez) == 0) return(NULL)
  
  uni_entrez <- NULL
  if (!is.null(universe)) {
    uni_entrez <- bitr(universe,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Hs.eg.db)$ENTREZID
  }
  
  enrichKEGG(
    gene          = entrez$ENTREZID,
    universe      = uni_entrez,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval,
    qvalueCutoff  = 0.2
  )
}

# ─────────────────────────────────────────────
# STEP 2: Run enrichment per contrast × tissue
#          separately for up and down genes
# ─────────────────────────────────────────────

# All background genes tested
universe_genes <- unique(de_overall_all$gene)

kegg_results <- de_overall_all %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) > 0.5,
    # Keep only clean comparisons — drop mirrored duplicates
    !contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010")
  ) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(tissue, contrast, direction) %>%
  group_map(~ {
    res <- run_kegg(.x$gene, universe = universe_genes)
    if (is.null(res)) return(NULL)
    
    as.data.frame(res) %>%
      mutate(
        tissue    = .y$tissue,
        contrast  = .y$contrast,
        direction = .y$direction
      )
  }, .keep = TRUE) %>%
  bind_rows()

# ─────────────────────────────────────────────
# STEP 3: Dot plot — top pathways per panel
# ─────────────────────────────────────────────

top_paths <- kegg_results %>%
  group_by(tissue, contrast, direction) %>%
  slice_min(p.adjust, n = 5) %>%
  ungroup()

ggplot(top_paths,
       aes(x     = direction,
           y     = reorder(Description, -p.adjust),
           size  = Count,
           color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#e63946", high = "#adb5bd",
                       name = "adj. p-value") +
  scale_size_continuous(range = c(2, 8), name = "Gene count") +
  facet_grid(cols = vars(tissue), rows = vars(contrast),
             scales = "free_y",   # ← each row gets its own y-axis
             space  = "free_y") + # ← panel height shrinks to fit its own genes
  theme_bw(base_size = 10) +
  theme(
    axis.text.y      = element_text(size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 8),
    legend.position  = "right"
  ) +
  labs(x = "Direction", y = NULL,
       title = "KEGG enrichment — vaccine comparisons by tissue")



# ─────────────────────────────────────────────
# STEP 3: Per-cluster DE — same design
# ─────────────────────────────────────────────
run_percluster_de <- function(seurat_obj,
                              cluster_col  = "RNA_snn_res.0.1",
                              contrast_var = "vaccine",
                              ref_level    = "Fluarix",
                              min_cells    = 10,
                              min_counts   = 10) {
  
  message("subsetting obj...")
  seurat_obj$pb_cluster <- seurat_obj[[cluster_col]][, 1]
  pb <- make_pseudobulk(
    seurat_obj,
    group_vars = c("patient_id", "tissue", "vaccine", "pb_cluster")
  )
  
  results_list <- list()
  
  for (tis in sort(unique(pb$meta$tissue))) {
    message("\n========================================")
    message("Processing tissue: ", tis)
    
    for (cl in sort(unique(pb$meta$pb_cluster))) {
      message("----------------------------------------")
      message("Cluster: ", cl)
      
      keep   <- pb$meta$tissue == tis & pb$meta$pb_cluster == cl
      counts <- pb$counts[, keep, drop = FALSE]
      s_meta <- pb$meta[keep, , drop = FALSE]
      
      message("Initial pseudobulk samples: ", ncol(counts))
      
      if (ncol(counts) < 4) {
        message("Skipping: < 4 samples"); next
      }
      
      keep2  <- colSums(counts) >= min_counts
      counts <- counts[, keep2, drop = FALSE]
      s_meta <- s_meta[keep2, , drop = FALSE]
      message("Samples after count filter (>= ", min_counts, "): ", ncol(counts))
      
      levels_present <- unique(s_meta[[contrast_var]])
      message("Levels present: ", paste(levels_present, collapse = ", "))
      
      # ── KEY FIX: skip if reference level not present ──
      if (!ref_level %in% levels_present) {
        message("Skipping: reference level '", ref_level, "' not present in this cluster/tissue")
        next
      }
      
      if (length(levels_present) < 2) {
        message("Skipping: only one level present"); next
      }
      
      s_meta[[contrast_var]] <- relevel(
        factor(s_meta[[contrast_var]]), ref = ref_level
      )
      
      keep_genes <- rowSums(counts >= 1) >= 2
      counts     <- counts[keep_genes, ]
      if (nrow(counts) < 10) {
        message("Skipping: too few genes"); next
      }
      
      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = counts,
          colData   = s_meta,
          design    = as.formula(paste("~", contrast_var))
        )
        dds <- DESeq(dds, quiet = TRUE)
        
        for (lvl in setdiff(levels(s_meta[[contrast_var]]), ref_level)) {
          res <- results(dds, contrast = c(contrast_var, lvl, ref_level)) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("gene") %>%
            mutate(tissue = tis, cluster = cl,
                   contrast = paste(lvl, "vs", ref_level))
          
          results_list[[paste(tis, cl, lvl, sep = "__")]] <- res
        }
        message("Finished tissue=", tis, ", cluster=", cl)
        
      }, error = function(e) {
        message("Error in tissue=", tis, " cluster=", cl, ": ", e$message)
      })
    }
  }
  
  bind_rows(results_list)
}
# Run it — will take a few minutes
de_percluster <- run_percluster_de(seurat_obj, ref_level = "Fluarix")
de_percluster_covid <- run_percluster_de(seurat_obj, ref_level = "covidmRNA")

de_percluster_all <- bind_rows(de_percluster, de_percluster_covid) %>%
  distinct(gene, tissue, cluster, contrast, .keep_all = TRUE) %>%  # drop mirrored duplicates
  filter(!contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010"))  # keep canonical direction only

summary_tbl <- de_percluster_all %>%
  filter(!is.na(padj)) %>%
  group_by(tissue, cluster, contrast) %>%
  summarise(
    n_sig_up   = sum(padj < 0.05 & log2FoldChange > 0.5),
    n_sig_down = sum(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  ) %>%
  arrange(tissue, cluster, contrast)

write.csv(
  summary_tbl,
  "de_percluster_summary.csv",
  row.names = FALSE
)

library(ggplot2)
library(dplyr)
library(ggrepel)

# Significant genes only, no mirrored contrasts
de_sig <- de_percluster_all %>%
  filter(
    !is.na(padj),
    !contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010")
  ) %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  0.5 ~ "up",
    padj < 0.05 & log2FoldChange < -0.5 ~ "down",
    TRUE                                 ~ "ns"
  ))

# Top genes to label per panel
top_genes <- de_sig %>%
  filter(sig != "ns") %>%
  group_by(tissue, contrast, cluster) %>%
  slice_min(padj, n = 10) %>%
  ungroup()

clusters <- sort(unique(de_sig$cluster))

cairo_pdf("DE_percluster_volcanos_by_contrast.pdf", width = 20, height = 16, onefile = TRUE)

contrasts <- unique(de_sig$contrast)
tissues   <- sort(unique(de_sig$tissue))
clusters  <- sort(unique(de_sig$cluster))

for (ct in contrasts) {
  for (tis in tissues) {
    
    dat      <- filter(de_sig,    contrast == ct, tissue == tis)
    dat_labs <- filter(top_genes, contrast == ct, tissue == tis)
    
    if (nrow(dat) == 0) next
    
    # Count sig genes per cluster panel
    counts <- dat %>%
      filter(sig != "ns") %>%
      count(cluster, sig) %>%
      tidyr::pivot_wider(names_from = sig, values_from = n, values_fill = 0) %>%
      mutate(
        up    = ifelse("up"   %in% colnames(.), up,   0),
        down  = ifelse("down" %in% colnames(.), down, 0),
        label = paste0("\u2191", up, "  \u2193", down)
      )
    
    p <- ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_text_repel(
        data         = dat_labs,
        aes(label    = gene),
        size         = 2.2,
        max.overlaps = 12,
        segment.size = 0.3,
        color        = "black"
      ) +
      geom_vline(xintercept = c(-0.5, 0.5),
                 linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_text(
        data = counts,
        aes(label = label, x = Inf, y = Inf),
        hjust = 1.1, vjust = 1.5,
        size = 3, color = "grey30",
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = c(up = "#e63946", down = "#457b9d", ns = "grey75"),
        labels = c(up = "Up", down = "Down", ns = "NS")
      ) +
      scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
      facet_wrap(~ cluster, ncol = 4,        # ← all clusters as panels
                 labeller = label_both) +
      theme_bw(base_size = 11) +
      theme(
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(size = 11, color = "grey40")
      ) +
      labs(
        x        = "log2 Fold Change",
        y        = "-log10(adjusted p-value)",
        color    = NULL,
        title    = paste0("Contrast: ", ct),
        subtitle = paste0("Tissue: ", tis, " — one panel per cluster")
      )
    
    print(p)
  }
}

dev.off()
message("Saved: DE_percluster_volcanos_by_contrast.pdf")

#############################################
# ─────────────────────────────────────────────
# MSigDB C2 enrichment helper
# ─────────────────────────────────────────────

library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)

msigdb_c2 <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
) %>%
  select(gs_name, gene_symbol)

run_msigdb_c2 <- function(genes, universe = NULL, pval = 0.05) {
  
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  
  if (length(genes) < 5) return(NULL)
  
  enr <- enricher(
    gene          = genes,
    universe      = universe,
    TERM2GENE     = msigdb_c2,
    pAdjustMethod = "BH",
    pvalueCutoff  = pval,
    qvalueCutoff  = 0.2
  )
  
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  
  return(enr)
}

# ─────────────────────────────────────────────
# Overall MSigDB C2 enrichment
# ─────────────────────────────────────────────

universe_genes_overall <- unique(de_overall_all$gene)

msigdb_c2_overall <- de_overall_all %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) > 0.5,
    !contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010")
  ) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(tissue, contrast, direction) %>%
  group_map(~ {
    
    res <- run_msigdb_c2(
      genes    = .x$gene,
      universe = universe_genes_overall
    )
    
    if (is.null(res)) return(NULL)
    
    as.data.frame(res) %>%
      mutate(
        tissue    = .y$tissue,
        contrast  = .y$contrast,
        direction = .y$direction
      )
    
  }, .keep = TRUE) %>%
  bind_rows()

write.csv(
  msigdb_c2_overall,
  "msigdb_c2_overall_enrichment.csv",
  row.names = FALSE
)

top_paths_overall <- msigdb_c2_overall %>%
  group_by(tissue, contrast, direction) %>%
  slice_min(p.adjust, n = 5) %>%
  ungroup()

ggplot(top_paths_overall,
       aes(x = direction,
           y = reorder(Description, -p.adjust),
           size = Count,
           color = p.adjust)) +
  geom_point() +
  scale_color_gradient(
    low = "#e63946",
    high = "#adb5bd",
    name = "adj. p-value"
  ) +
  scale_size_continuous(range = c(2, 8), name = "Gene count") +
  facet_grid(
    cols = vars(tissue),
    rows = vars(contrast),
    scales = "free_y",
    space = "free_y"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y      = element_text(size = 8),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 8),
    legend.position  = "right"
  ) +
  labs(
    x = "Direction",
    y = NULL,
    title = "MSigDB C2 enrichment — overall vaccine comparisons by tissue"
  )
merged_after_reclustering.rds


# Save results
saveRDS(de_overall,    "pseudobulk_DE_overall.rds")
saveRDS(de_percluster, "pseudobulk_DE_percluster.rds")
