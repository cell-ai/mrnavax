###crid

library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)
library(Seurat)

setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

# ── paths ─────────────────────────────────────────────────────────
recluster_dir <- "./markers_0.1_reclustered"
bt_dir        <- "./markers_BT_multiresolution"
plot_dir      <- "./plots_final_annotations"
dir.create(plot_dir, showWarnings = FALSE)

# ── color palettes ────────────────────────────────────────────────
palette_myeloid <- c(
  "Classical CD14+ monocyte"        = "#E63946",
  "pDC"                             = "#457B9D",
  "cDC2 resident"                   = "#2A9D8F",
  "Non-classical CD16+ monocyte"    = "#E9C46A",
  "Migratory cDC2 CCR7+"            = "#F4A261",
  "T cell contamination"            = "#CCCCCC",
  "Mast cell progenitor"            = "#8338EC",
  "Pro-B lymphoid progenitor"       = "#3A86FF"
)

myeloid_map <- c(
  "0" = "Classical CD14+ monocyte",
  "1" = "pDC",
  "2" = "cDC2 resident",
  "3" = "Non-classical CD16+ monocyte",
  "4" = "Migratory cDC2 CCR7+",
  "5" = "T cell contamination",
  "6" = "Mast cell progenitor",
  "7" = "Pro-B lymphoid progenitor"
)

nk_map <- c(
  "0" = "T / NKT cell",
  "1" = "Cytotoxic NK cell",
  "2" = "Quiescent NK / ILC"
)

palette_nk <- c(
  "T / NKT cell"       = "#CCCCCC",
  "Cytotoxic NK cell"  = "#E63946",
  "Quiescent NK / ILC" = "#457B9D"
)

# plasma map is correct — cluster 3 has only 19 cells
plasma_map <- c(
  "0" = "Plasmablast",
  "1" = "Long-lived plasma cell",
  "2" = "GC B cell",
  "3" = "T cell contamination"
)

palette_plasma <- c(
  "Plasmablast"            = "#E63946",
  "Long-lived plasma cell" = "#F4A261",
  "GC B cell"              = "#2A9D8F",
  "T cell contamination"   = "#CCCCCC"
)

# B map confirmed correct against cluster numbers
b_map <- c(
  "0" = "Non-B contamination",
  "1" = "IgA plasma cell",
  "2" = "Transitional B cell",
  "3" = "Early activated B cell",
  "4" = "IgA memory B cell",
  "5" = "Naive B cell",
  "6" = "Pre-B progenitor"
)

palette_b <- c(
  "Non-B contamination"    = "#CCCCCC",
  "IgA plasma cell"        = "#E63946",
  "Transitional B cell"    = "#F4A261",
  "Early activated B cell" = "#2A9D8F",
  "IgA memory B cell"      = "#8338EC",
  "Naive B cell"           = "#457B9D",
  "Pre-B progenitor"       = "#3A86FF",
  "IFN-stimulated B cell"  = "hotpink3"
)

palette_t <- c(
  "CD8+ effector T cell"              = "#E63946",
  "Early activated T cell"            = "#F4A261",
  "CD4+ effector memory T cell"       = "#2A9D8F",
  "Naive T cell"                      = "#457B9D",
  "IFN-stimulated T cell"             = "#FB8500",
  "Quiescent / central memory T cell" = "#8338EC",
  "Proliferating T cell"              = "#3A86FF"
)

# ── helper: extract umap + metadata from SCE ─────────────────────
sce_to_df <- function(sce, umap_key = "X_umap") {
  
  # get UMAP — zellkonverter stores obsm as reducedDims
  rd_names <- reducedDimNames(sce)
  message("  available reducedDims: ", paste(rd_names, collapse = ", "))
  
  # try common names
  umap_name <- intersect(
    c("X_umap", "UMAP", "umap", "X_UMAP", "UMAP_SUB"),
    rd_names
  )[1]
  
  if (is.na(umap_name)) {
    stop("No UMAP found in reducedDims. Available: ", paste(rd_names, collapse = ", "))
  }
  
  message("  using reducedDim: ", umap_name)
  umap <- as.data.frame(reducedDim(sce, umap_name))
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  
  meta <- as.data.frame(colData(sce))
  df   <- cbind(meta, umap)
  return(df)
}

# ── helper: UMAP plot ─────────────────────────────────────────────
plot_umap_df <- function(df, color_col, title, palette,
                         pt_size = 0.05, label = TRUE) {
  
  labels <- df %>%
    group_by(.data[[color_col]]) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2,
                      color = .data[[color_col]])) +
    geom_point(size = pt_size, alpha = 0.4, stroke = 0) +
    scale_color_manual(values = palette, na.value = "#DDDDDD") +
    ggtitle(title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.title = element_blank(),
      legend.text  = element_text(size = 9),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.line    = element_blank(),
      axis.title   = element_text(size = 9, color = "grey50"),
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 4, alpha = 1)
    ))
  
  if (label) {
    p <- p + geom_label_repel(
      data         = labels,
      aes(label    = .data[[color_col]]),
      color        = "black",
      fill         = "white",
      size         = 2.8,
      alpha        = 0.85,
      label.size   = 0.2,
      max.overlaps = 20,
      show.legend  = FALSE
    )
  }
  
  return(p)
}

# ══════════════════════════════════════════════════════════════════
# MYELOID
# ══════════════════════════════════════════════════════════════════
message("\n── Myeloid ──")



sce <- readH5AD(file.path(recluster_dir, "Myeloid_scgpt.h5ad"))
df  <- sce_to_df(sce)
# check what values exist in the reclustered column
table(colData(sce)$leiden_res_0.1_reclustered)

# update the map accordingly
df$cell_type <- myeloid_map[as.character(df$leiden_res_0.1_reclustered)]

# and in the plot function call
df_clean <- df %>%
  filter(!is.na(cell_type),
         cell_type != "T cell contamination")

p_myeloid <- plot_umap_df(
  df_clean, "cell_type",
  "Myeloid cells", palette_myeloid
)
# ══════════════════════════════════════════════════════════════════
# NK
# ══════════════════════════════════════════════════════════════════
message("\n── NK ──")

sce    <- readH5AD(file.path(recluster_dir, "NK_scgpt.h5ad"))
df     <- sce_to_df(sce)

df$cell_type <- nk_map[as.character(df$leiden_res_0.1_reclustered)]

df_clean <- df %>% filter(!is.na(cell_type),
                          cell_type != "T cell contamination")
print(df_clean %>% count(cell_type) %>% arrange(desc(n)))

p_nk <- plot_umap_df(df_clean, "cell_type", "NK cells", palette_nk)
ggsave(file.path(plot_dir, "umap_NK_annotated.pdf"),
       p_nk, width = 9, height = 6)
message("saved umap_NK_annotated.pdf")
rm(sce, df, df_clean); gc()

# ══════════════════════════════════════════════════════════════════
# PLASMA
# ══════════════════════════════════════════════════════════════════
message("\n── Plasma ──")

sce <- readH5AD(file.path(recluster_dir, "Plasma_scgpt.h5ad"))
df  <- sce_to_df(sce)
# and make sure you use the reclustered column
df$cell_type <- plasma_map[as.character(df$leiden_res_0.1_reclustered)]

df_clean <- df %>% filter(!is.na(cell_type))

print(df_clean %>% count(cell_type) %>% arrange(desc(n)))

p_plasma <- plot_umap_df(df_clean, "cell_type", "Plasma cells", palette_plasma)
ggsave(file.path(plot_dir, "umap_Plasma_annotated.pdf"),
       p_plasma, width = 9, height = 6)
rm(sce, df); gc()

# ══════════════════════════════════════════════════════════════════
# B CELLS
# ══════════════════════════════════════════════════════════════════
message("\n── B cells ──")

sce <- readH5AD(file.path(bt_dir, "B_scgpt.h5ad"))
df  <- sce_to_df(sce)

df$cell_type <- b_map[as.character(df$leiden_res_0.3)]
df$cell_type <- ifelse(
  as.character(df$leiden_res_0.5) == "4",
  "IFN-stimulated B cell",
  df$cell_type
)

p_b <- plot_umap_df(df, "cell_type", "B cells", palette_b)
ggsave(file.path(plot_dir, "umap_B_annotated.pdf"),
       p_b, width = 10, height = 7)
message("saved umap_B_annotated.pdf")
rm(sce, df); gc()

library(UCell)

# compute IFN module score on B cells
ifn_genes <- c("IFITM1", "IFITM2", "MX1", "STAT1", 
               "ISG15", "ISG20", "LY6E", "IFI6", "TRIM22")

# add as continuous score overlaid on UMAP
# first get the SCE object back
sce_b <- readH5AD(file.path(bt_dir, "B_scgpt.h5ad"))

scores <- ScoreSignatures_UCell(
  assay(sce_b, "logcounts"),
  features = list(IFN_response = ifn_genes)
)

df$IFN_score <- as.numeric(scores)

# plot as continuous overlay
p_b_ifn <- ggplot(df %>% arrange(IFN_score),
                  aes(x = UMAP_1, y = UMAP_2, color = IFN_score)) +
  geom_point(size = 0.05, alpha = 0.5, stroke = 0) +
  scale_color_gradientn(
    colors = c("lightgrey", "#FED976", "#FD8D3C", "#E31A1C"),
    name   = "IFN score"
  ) +
  ggtitle("B cells — IFN response score") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank()
  )

ggsave(file.path(plot_dir, "umap_B_IFN_score.pdf"),
       p_b_ifn, width = 8, height = 6)
# ══════════════════════════════════════════════════════════════════
# T CELLS
# ══════════════════════════════════════════════════════════════════
message("\n── T cells ──")

sce <- readH5AD(file.path(bt_dir, "T_scgpt_annotated.h5ad"))
df  <- sce_to_df(sce)
df$cell_type <- df$cell_type_T

p_t <- plot_umap_df(df, "cell_type", "T cells", palette_t, pt_size = 0.03)
ggsave(file.path(plot_dir, "umap_T_annotated.pdf"),
       p_t, width = 11, height = 7)
message("saved umap_T_annotated.pdf")
rm(sce, df); gc()

library(UCell)

naive_genes    <- c("LEF1", "SATB1", "CCR7", "SELL", "TCF7", "IL7R", "NELL2")
effector_genes <- c("GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG", "TNF")
ifn_genes      <- c("MX1", "STAT1", "ISG15", "IFITM1", "IFI6", "IFI44L", "OAS1")

sce_t <- readH5AD(file.path(bt_dir, "T_scgpt_annotated.h5ad"))

scores <- ScoreSignatures_UCell(
  assay(sce_t, "logcounts"),
  features = list(
    Naive_score    = naive_genes,
    Effector_score = effector_genes,
    IFN_score      = ifn_genes
  )
)

df_t <- as.data.frame(colData(sce_t))
df_t$Naive_score    <- scores[, "Naive_score_UCell"]
df_t$Effector_score <- scores[, "Effector_score_UCell"]
df_t$IFN_score      <- scores[, "IFN_score_UCell"]

# compare scores within IFN-stimulated vs other clusters
df_t %>%
  group_by(cell_type_T) %>%
  summarise(
    mean_naive    = mean(Naive_score),
    mean_effector = mean(Effector_score),
    mean_ifn      = mean(IFN_score)
  ) %>%
  arrange(desc(mean_ifn))

df_t %>%
  filter(!is.na(cell_type_T)) %>%
  group_by(timepoint, cell_type_T) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  filter(cell_type_T == "IFN-stimulated T cell") %>%
  arrange(timepoint) %>%
  print(n = Inf)

ifn_props <- df_t %>%
  mutate(lineage = "T") %>%
  group_by(donor, timepoint, vaccine, lineage) %>%
  summarise(
    ifn_prop = mean(cell_type_T == "IFN-stimulated T cell"),
    .groups = "drop"
  )

# check structure
print(ifn_props %>% arrange(vaccine, timepoint))

# proportion of IFN-stimulated T cells per timepoint per vaccine
ifn_props %>%
  group_by(vaccine, timepoint) %>%
  summarise(
    mean_prop = mean(ifn_prop),
    sd_prop   = sd(ifn_prop),
    n_donors  = n(),
    .groups   = "drop"
  ) %>%
  arrange(vaccine, timepoint) %>%
  print(n = Inf)

library(ggplot2)

# convert timepoint to numeric for plotting
timepoint_days <- c(
  "d0" = 0, "d8" = 8, "d15" = 15, "d28" = 28,
  "d28+d35" = 31, "d35" = 35, "d57" = 57, "d60" = 60,
  "d110" = 110, "d121" = 121, "d180" = 180,
  "d181" = 181, "d201" = 201
)

ifn_summary <- ifn_props %>%
  group_by(vaccine, timepoint) %>%
  summarise(
    mean_prop = mean(ifn_prop),
    sd_prop   = sd(ifn_prop),
    n_donors  = n(),
    se_prop   = sd_prop / sqrt(n_donors),
    .groups   = "drop"
  ) %>%
  mutate(day = timepoint_days[as.character(timepoint)])

ggplot(ifn_summary, aes(x = day, y = mean_prop * 100,
                        color = vaccine, group = vaccine)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_donors)) +
  geom_errorbar(
    aes(ymin = (mean_prop - se_prop) * 100,
        ymax = (mean_prop + se_prop) * 100),
    width = 3
  ) +
  scale_color_manual(values = c(
    "Fluarix"   = "#E63946",
    "covidmRNA" = "#457B9D",
    "mRNA-1010" = "#2A9D8F"
  )) +
  scale_size_continuous(range = c(2, 5), name = "n donors") +
  labs(
    x     = "Days post vaccination",
    y     = "IFN-stimulated T cells (%)",
    color = "Vaccine",
    title = "IFN-stimulated T cell kinetics post vaccination"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5)

ggsave(file.path(plot_dir, "IFN_T_kinetics_by_vaccine.pdf"),
       width = 10, height = 6)
# ══════════════════════════════════════════════════════════════════
# COMBINED PANEL
# ══════════════════════════════════════════════════════════════════
message("\n── combined panel ──")

combined <- (p_myeloid + p_nk)  /
  (p_plasma  + p_b)   /
  (p_t       + plot_spacer()) +
  plot_annotation(
    title = "Single-cell annotation — all lineages",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
    )
  )

ggsave(
  file.path(plot_dir, "umap_all_lineages_panel.pdf"),
  combined, width = 20, height = 24
)
message("saved umap_all_lineages_panel.pdf")

message("\n── all done ──")
