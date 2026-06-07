setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

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

select <- dplyr::select
filter <- dplyr::filter

# ── 1. Load object ─────────────────────────────────────────────────────────────
seurat_obj <- readRDS("merged_after_reclustering_no_TCR_BCR.rds")

cluster_to_celltype <- c(
  "0"  = "Naive/central memory CD4 T",
  "1"  = "Activated/memory B cells",
  "2"  = "Naive B",
  "3"  = "Plasma cells",
  "4"  = "Cytotoxic CD8 T",
  "5"  = "Germinal center B",
  "6"  = "Cytotoxic lymphocytes",
  "7"  = "T follicular helper (Tfh)",
  "8"  = "Classical monocytes",
  "9"  = "Activated/memory T",
  "10" = "Plasmablasts",
  "11" = "Plasmacytoid dendritic cells (pDC)"
)

cell_type <- cluster_to_celltype[as.character(seurat_obj$RNA_snn_res.0.1)]
names(cell_type) <- colnames(seurat_obj)
seurat_obj <- AddMetaData(seurat_obj, metadata = cell_type, col.name = "cell_type")

# ── 2. Assign timepoint phase ──────────────────────────────────────────────────
early_pbmc <- c("d0", "d8", "d28", "d28+d35")
late_pbmc  <- c("d121", "d180")

early_ln   <- c("d0", "d15", "d28", "d35")
late_ln    <- c("d57", "d60", "d110", "d121", "d181", "d201")

seurat_obj$timepoint_phase <- case_when(
  seurat_obj$tissue == "PBMC"       &
    seurat_obj$timepoint %in% early_pbmc ~ "early",
  seurat_obj$tissue == "PBMC"       &
    seurat_obj$timepoint %in% late_pbmc  ~ "late",
  seurat_obj$tissue == "lymph_node" &
    seurat_obj$timepoint %in% early_ln   ~ "early",
  seurat_obj$tissue == "lymph_node" &
    seurat_obj$timepoint %in% late_ln    ~ "late",
  TRUE ~ NA_character_
)

# Verify
cat("Timepoint phase distribution:\n")
print(table(seurat_obj$tissue,
            seurat_obj$timepoint_phase,
            useNA = "ifany"))

cat("\nPhase × vaccine × tissue:\n")
print(
  seurat_obj@meta.data %>%
    dplyr::distinct(patient_id, tissue, vaccine,
                    timepoint, timepoint_phase) %>%
    dplyr::count(tissue, vaccine, timepoint_phase, timepoint) %>%
    dplyr::arrange(tissue, vaccine, timepoint_phase, timepoint)
)

# ── 3. MSigDB C2 setup ─────────────────────────────────────────────────────────
get_msigdb_c2 <- function() {
  args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% args) {
    msigdbr(species = "Homo sapiens", collection = "C2")
  } else {
    msigdbr(species = "Homo sapiens", category = "C2")
  }
}

msigdb_c2 <- get_msigdb_c2() %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

# ── 4. Helpers ─────────────────────────────────────────────────────────────────
clean_contrasts <- function(de_tbl) {
  de_tbl %>%
    distinct(gene, tissue, across(any_of(c("cluster", "timepoint_phase"))),
             contrast, .keep_all = TRUE) %>%
    filter(!contrast %in% c("Fluarix vs covidmRNA", "Fluarix vs mRNA-1010"))
}

make_pseudobulk <- function(obj, group_vars, assay = "RNA") {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")
  
  meta <- obj[[]] %>%
    filter(!is.na(patient_id), !is.na(vaccine))
  
  valid_cells <- rownames(meta)[
    complete.cases(meta[, group_vars, drop = FALSE])
  ]
  meta   <- meta[valid_cells, , drop = FALSE]
  counts <- counts[, valid_cells, drop = FALSE]
  
  meta$pb_group <- apply(
    meta[, group_vars, drop = FALSE], 1, paste, collapse = "__"
  )
  
  groups <- unique(meta$pb_group)
  
  pb_counts <- sapply(groups, function(g) {
    cells <- rownames(meta)[meta$pb_group == g]
    if (length(cells) == 1) as.numeric(counts[, cells])
    else Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  rownames(pb_counts) <- rownames(counts)
  colnames(pb_counts) <- groups
  
  pb_meta <- meta %>%
    rownames_to_column("cell") %>%
    group_by(pb_group) %>%
    slice(1) %>%
    ungroup() %>%
    column_to_rownames("pb_group") %>%
    dplyr::select(all_of(group_vars))
  
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  list(counts = pb_counts, meta = pb_meta)
}

run_deseq <- function(counts, s_meta,
                      contrast_var, ref_level,
                      extra_info,
                      min_counts = 10) {
  
  keep_samples <- colSums(counts) >= min_counts
  counts  <- counts[, keep_samples, drop = FALSE]
  s_meta  <- s_meta[keep_samples, , drop = FALSE]
  
  if (ncol(counts) < 4)               return(NULL)
  levels_present <- unique(s_meta[[contrast_var]])
  if (!ref_level %in% levels_present) return(NULL)
  if (length(levels_present) < 2)     return(NULL)
  
  s_meta[[contrast_var]] <- relevel(
    factor(s_meta[[contrast_var]]), ref = ref_level
  )
  
  keep_genes <- rowSums(counts >= 1) >= 2
  counts <- counts[keep_genes, , drop = FALSE]
  if (nrow(counts) < 10) return(NULL)
  
  tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = round(counts),
      colData   = s_meta,
      design    = as.formula(paste("~", contrast_var))
    )
    dds <- DESeq(dds, quiet = TRUE)
    
    bind_rows(lapply(
      setdiff(levels(s_meta[[contrast_var]]), ref_level),
      function(lvl) {
        results(dds, contrast = c(contrast_var, lvl, ref_level)) %>%
          as.data.frame() %>%
          rownames_to_column("gene") %>%
          mutate(contrast = paste(lvl, "vs", ref_level), !!!extra_info)
      }
    ))
  }, error = function(e) {
    message("  DESeq2 error: ", e$message); NULL
  })
}

run_msigdb_c2_enrich <- function(genes, universe = NULL, pval = 0.05) {
  genes <- unique(genes[!is.na(genes)])
  if (length(genes) < 5) return(NULL)
  if (!is.null(universe)) universe <- unique(universe[!is.na(universe)])
  
  enr <- enricher(
    gene          = genes,
    universe      = universe,
    TERM2GENE     = msigdb_c2,
    pAdjustMethod = "BH",
    pvalueCutoff  = pval,
    qvalueCutoff  = 0.2
  )
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  enr
}

run_kegg_enrich <- function(genes, universe = NULL, pval = 0.05) {
  genes <- unique(genes[!is.na(genes)])
  if (length(genes) < 5) return(NULL)
  
  entrez <- suppressMessages(
    bitr(genes, fromType = "SYMBOL",
         toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(entrez) || nrow(entrez) < 5) return(NULL)
  
  uni_entrez <- NULL
  if (!is.null(universe)) {
    uni_entrez <- suppressMessages(
      bitr(unique(universe), fromType = "SYMBOL",
           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )$ENTREZID
  }
  
  enr <- enrichKEGG(
    gene          = unique(entrez$ENTREZID),
    universe      = unique(uni_entrez),
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval,
    qvalueCutoff  = 0.2
  )
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  enr
}

# Generic enrichment runner
run_enrichment <- function(de_tbl, group_cols,
                           output_prefix,
                           universe    = NULL,
                           padj_cutoff = 0.05,
                           lfc_cutoff  = 0.5) {
  
  de_sig <- de_tbl %>%
    filter(!is.na(padj), !is.na(log2FoldChange),
           padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down"))
  
  if (nrow(de_sig) == 0) {
    message("  No significant genes for: ", output_prefix)
    return(list(msigdb = NULL, kegg = NULL))
  }
  
  run_one <- function(run_fun, label) {
    results <- de_sig %>%
      group_by(across(all_of(c(group_cols, "direction")))) %>%
      group_map(~ {
        res <- run_fun(.x$gene, universe = universe)
        if (is.null(res)) return(NULL)
        as.data.frame(res) %>%
          mutate(!!!as.list(.y),
                 gene_input_n = n_distinct(.x$gene))
      }, .keep = TRUE) %>%
      bind_rows()
    
    if (nrow(results) == 0) return(NULL)
    
    fname <- paste0(output_prefix, "_", label)
    write.csv(results,  paste0(fname, ".csv"),  row.names = FALSE)
    saveRDS(results,    paste0(fname, ".rds"))
    message("  Saved: ", fname)
    results
  }
  
  list(
    msigdb = run_one(run_msigdb_c2_enrich, "msigdb_c2"),
    kegg   = run_one(run_kegg_enrich,      "kegg")
  )
}

# ── 5. Volcano plot helper ─────────────────────────────────────────────────────
save_volcano <- function(de_tbl, facet_col, filename,
                         title = "", width = 20, height = 14) {
  
  de_plot <- de_tbl %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(sig = case_when(
      padj < 0.05 & log2FoldChange >  0.5 ~ "up",
      padj < 0.05 & log2FoldChange < -0.5 ~ "down",
      TRUE ~ "ns"
    ))
  
  top_labs <- de_plot %>%
    filter(sig != "ns") %>%
    group_by(across(all_of(c(facet_col, "contrast", "tissue")))) %>%
    slice_min(padj, n = 8) %>%
    ungroup()
  
  p <- ggplot(de_plot,
              aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_text_repel(data = top_labs, aes(label = gene),
                    size = 2, max.overlaps = 12,
                    segment.size = 0.3, colour = "black") +
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(
      values = c(up = "#e63946", down = "#457b9d", ns = "grey75")) +
    scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
    facet_grid(rows = vars(contrast),
               cols = vars(.data[[facet_col]]),
               scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92"),
          strip.text        = element_text(face = "bold", size = 8),
          legend.position   = "bottom",
          panel.grid.minor  = element_blank()) +
    labs(x = "log2FC", y = "-log10(adj p)", colour = NULL,
         title = title)
  
  ggsave(filename, p, width = width, height = height)
  message("Saved: ", filename)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. OVERALL DE — stratified by tissue × timepoint_phase × vaccine
# ══════════════════════════════════════════════════════════════════════════════

message("\n════ OVERALL DE ════")

run_overall_de_phase <- function(obj, ref_level,
                                 contrast_var = "vaccine") {
  pb <- make_pseudobulk(
    obj,
    group_vars = c("patient_id", "tissue", "vaccine", "timepoint_phase")
  )
  
  bind_rows(lapply(
    sort(unique(pb$meta$tissue)), function(tis) {
      bind_rows(lapply(
        sort(unique(pb$meta$timepoint_phase)), function(phase) {
          if (is.na(phase)) return(NULL)
          
          message("  Overall DE | tissue=", tis,
                  " | phase=", phase,
                  " | ref=", ref_level)
          
          keep   <- pb$meta$tissue == tis &
            !is.na(pb$meta$timepoint_phase) &
            pb$meta$timepoint_phase == phase
          counts <- pb$counts[, keep, drop = FALSE]
          s_meta <- pb$meta[keep,  , drop = FALSE]
          
          run_deseq(counts, s_meta,
                    contrast_var = contrast_var,
                    ref_level    = ref_level,
                    extra_info   = list(tissue          = tis,
                                        timepoint_phase = phase))
        })
      ))
    }
  ))
}

de_overall_fluarix <- run_overall_de_phase(seurat_obj, "Fluarix")
de_overall_covid   <- run_overall_de_phase(seurat_obj, "covidmRNA")

de_overall_all <- bind_rows(de_overall_fluarix, de_overall_covid) %>%
  clean_contrasts()

write.csv(de_overall_all, "enrichment/de_overall_phase_all_contrasts.csv",
          row.names = FALSE)

# Summary
overall_summary <- de_overall_all %>%
  filter(!is.na(padj)) %>%
  group_by(tissue, timepoint_phase, contrast) %>%
  summarise(
    n_up   = sum(padj < 0.05 & log2FoldChange >  0.5),
    n_down = sum(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  )
write.csv(overall_summary, "enrichment/de_overall_phase_summary.csv", row.names = FALSE)
print(overall_summary)

# Volcano — one column per tissue × phase
cairo_pdf("volcano_overall_by_phase.pdf",
          width = 20, height = 10, onefile = TRUE)
for (phase in c("early", "late")) {
  dat <- de_overall_all %>%
    filter(timepoint_phase == phase,
           !is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(sig = case_when(
      padj < 0.05 & log2FoldChange >  0.5 ~ "up",
      padj < 0.05 & log2FoldChange < -0.5 ~ "down",
      TRUE ~ "ns"
    ))
  if (nrow(dat) == 0) next
  
  top_labs <- dat %>%
    filter(sig != "ns") %>%
    group_by(tissue, contrast) %>%
    slice_min(padj, n = 8) %>%
    ungroup()
  
  p <- ggplot(dat, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_text_repel(data = top_labs, aes(label = gene),
                    size = 2, max.overlaps = 12,
                    segment.size = 0.3, colour = "black") +
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(
      values = c(up = "#e63946", down = "#457b9d", ns = "grey75")) +
    scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
    facet_grid(rows = vars(contrast), cols = vars(tissue)) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92"),
          strip.text        = element_text(face = "bold", size = 9),
          legend.position   = "bottom") +
    labs(title    = paste0("Overall DE — ", phase, " phase"),
         x = "log2FC", y = "-log10(adj p)", colour = NULL)
  print(p)
}
dev.off()

# Overall enrichment — stratified by tissue × phase × contrast
message("\nOverall enrichment (MSigDB + KEGG)...")

enrich_overall <- run_enrichment(
  de_tbl        = de_overall_all,
  group_cols    = c("tissue", "timepoint_phase", "contrast"),
  output_prefix = "enrich_overall_phase",
  universe      = unique(de_overall_all$gene)
)

# ══════════════════════════════════════════════════════════════════════════════
# 7. PER-CLUSTER DE — stratified by tissue × timepoint_phase × vaccine
# ══════════════════════════════════════════════════════════════════════════════

message("\n════ PER-CLUSTER DE ════")

run_percluster_de_phase <- function(obj,
                                    ref_level,
                                    cluster_col  = "RNA_snn_res.0.1",
                                    contrast_var = "vaccine") {
  obj$pb_cluster <- obj[[cluster_col]][, 1]
  
  pb <- make_pseudobulk(
    obj,
    group_vars = c("patient_id", "tissue", "vaccine",
                   "timepoint_phase", "pb_cluster")
  )
  
  bind_rows(lapply(
    sort(unique(pb$meta$tissue)), function(tis) {
      bind_rows(lapply(
        sort(unique(pb$meta$timepoint_phase)), function(phase) {
          if (is.na(phase)) return(NULL)
          bind_rows(lapply(
            sort(unique(pb$meta$pb_cluster)), function(cl) {
              
              message("  Per-cluster DE | tissue=", tis,
                      " | phase=", phase,
                      " | cluster=", cl,
                      " | ref=", ref_level)
              
              keep <- pb$meta$tissue == tis &
                !is.na(pb$meta$timepoint_phase) &
                pb$meta$timepoint_phase == phase &
                pb$meta$pb_cluster == cl
              
              counts <- pb$counts[, keep, drop = FALSE]
              s_meta <- pb$meta[keep,  , drop = FALSE]
              
              run_deseq(counts, s_meta,
                        contrast_var = contrast_var,
                        ref_level    = ref_level,
                        extra_info   = list(
                          tissue          = tis,
                          timepoint_phase = phase,
                          cluster         = cl,
                          cell_type       = unname(
                            cluster_to_celltype[as.character(cl)]
                          )
                        ))
            })
          ))
        })
  ))
    }
  ))
}

de_pc_fluarix <- run_percluster_de_phase(seurat_obj, "Fluarix")
de_pc_covid   <- run_percluster_de_phase(seurat_obj, "covidmRNA")

de_pc_all <- bind_rows(de_pc_fluarix, de_pc_covid) %>%
  clean_contrasts()

write.csv(de_pc_all, "enrichment/de_percluster_phase_all_contrasts.csv",
          row.names = FALSE)

# Summary
pc_summary <- de_pc_all %>%
  filter(!is.na(padj)) %>%
  group_by(tissue, timepoint_phase, cluster, cell_type, contrast) %>%
  summarise(
    n_up   = sum(padj < 0.05 & log2FoldChange >  0.5),
    n_down = sum(padj < 0.05 & log2FoldChange < -0.5),
    .groups = "drop"
  ) %>%
  arrange(tissue, timepoint_phase, cluster, contrast)

write.csv(pc_summary, "enrichment/de_percluster_phase_summary.csv", row.names = FALSE)
print(pc_summary)

# Per-cluster volcano PDFs — one page per phase × tissue
cairo_pdf("volcano_percluster_by_phase.pdf",
          width = 22, height = 16, onefile = TRUE)

for (phase in c("early", "late")) {
  for (tis in sort(unique(de_pc_all$tissue))) {
    for (ct in unique(de_pc_all$contrast)) {
      
      dat <- de_pc_all %>%
        filter(timepoint_phase == phase,
               tissue == tis, contrast == ct,
               !is.na(padj), !is.na(log2FoldChange)) %>%
        mutate(sig = case_when(
          padj < 0.05 & log2FoldChange >  0.5 ~ "up",
          padj < 0.05 & log2FoldChange < -0.5 ~ "down",
          TRUE ~ "ns"
        ))
      
      if (nrow(dat) == 0) next
      
      top_labs <- dat %>%
        filter(sig != "ns") %>%
        group_by(cluster) %>%
        slice_min(padj, n = 8) %>%
        ungroup()
      
      counts_lab <- dat %>%
        filter(sig != "ns") %>%
        count(cluster, sig) %>%
        pivot_wider(names_from = sig, values_from = n,
                    values_fill = 0) %>%
        mutate(
          up   = if ("up"   %in% colnames(.)) up   else 0L,
          down = if ("down" %in% colnames(.)) down else 0L,
          label = paste0("\u2191", up, "  \u2193", down)
        )
      
      p <- ggplot(dat,
                  aes(log2FoldChange, -log10(padj), colour = sig)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_text_repel(data = top_labs, aes(label = gene),
                        size = 2, max.overlaps = 12,
                        segment.size = 0.3, colour = "black") +
        geom_vline(xintercept = c(-0.5, 0.5),
                   linetype = "dashed", colour = "grey50",
                   linewidth = 0.3) +
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed", colour = "grey50",
                   linewidth = 0.3) +
        geom_text(data = counts_lab,
                  aes(label = label, x = Inf, y = Inf),
                  hjust = 1.1, vjust = 1.5,
                  size = 3, colour = "grey30",
                  inherit.aes = FALSE) +
        scale_colour_manual(
          values = c(up = "#e63946", down = "#457b9d", ns = "grey75")) +
        scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
        facet_wrap(~ cluster + cell_type, ncol = 4,
                   labeller = label_both) +
        theme_bw(base_size = 10) +
        theme(strip.background = element_rect(fill = "grey92"),
              strip.text        = element_text(face = "bold", size = 7),
              legend.position   = "bottom",
              panel.grid.minor  = element_blank()) +
        labs(
          title    = paste0("Per-cluster DE | ",
                            phase, " | ", ct),
          subtitle = paste0("Tissue: ", tis),
          x = "log2FC", y = "-log10(adj p)", colour = NULL
        )
      print(p)
    }
  }
}
dev.off()
message("Saved: volcano_percluster_by_phase.pdf")

# Per-cluster enrichment
message("\nPer-cluster enrichment (MSigDB + KEGG)...")

enrich_percluster <- run_enrichment(
  de_tbl        = de_pc_all,
  group_cols    = c("tissue", "timepoint_phase", "cluster", "contrast"),
  output_prefix = "enrich_percluster_phase",
  universe      = unique(de_pc_all$gene)
)

# ── 8. Save ────────────────────────────────────────────────────────────────────
saveRDS(de_overall_all,    "enrichment/de_overall_phase.rds")
saveRDS(de_pc_all,         "enrichment/de_percluster_phase.rds")
saveRDS(enrich_overall,    "enrichment/enrich_overall_phase.rds")
saveRDS(enrich_percluster, "enrichment/enrich_percluster_phase.rds")

message("\nDone.")
