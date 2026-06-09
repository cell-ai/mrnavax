# =============================================================================
# Lineage Cluster Analysis
# Analyses:
#   1. Clusters unique to one vaccine (100% of cells from a single vaccine)
#   2. Clusters appearing in only one patient (patient_id)
#   3. Cell count distribution per cluster within each lineage
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

RDS_PATH    <- "results/subobjects_by_lineage.rds"
OUTPUT_DIR  <- "results/lineage_analysis"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Resolution to use per lineage
res_map <- c(
  T       = "RNA_snn_res.0.3",
  NK      = "RNA_snn_res.0.1",
  B       = "RNA_snn_res.0.1",
  Plasma  = "RNA_snn_res.0.05",
  Myeloid = "RNA_snn_res.0.1"
)

message("Loading RDS...")
obj_list <- readRDS(RDS_PATH)
message("Loaded lineages: ", paste(names(obj_list), collapse = ", "))

# -----------------------------------------------------------------------------
# Helper: extract metadata with chosen resolution as 'cluster'
# -----------------------------------------------------------------------------

get_meta <- function(sobj, res_col) {
  md <- sobj@meta.data
  md$cluster <- as.character(md[[res_col]])
  md$cell_barcode <- rownames(md)
  md
}

# =============================================================================
# Analysis 1 — Clusters 100% unique to one vaccine
# =============================================================================

message("\n--- Analysis 1: Vaccine-unique clusters ---")

vaccine_unique_list <- lapply(names(obj_list), function(lin) {
  md  <- get_meta(obj_list[[lin]], res_map[lin])
  
  tbl <- md %>%
    group_by(cluster, vaccine) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(
      total_cells   = sum(n_cells),
      n_vaccines    = n_distinct(vaccine),
      pct           = n_cells / total_cells * 100
    ) %>%
    ungroup()
  
  # Keep only clusters where ALL cells belong to exactly one vaccine
  unique_clusters <- tbl %>%
    filter(n_vaccines == 1) %>%
    dplyr::select(cluster, vaccine, n_cells, total_cells, pct) %>%
    mutate(lineage = lin)
  
  unique_clusters
})

vaccine_unique_df <- bind_rows(vaccine_unique_list)

message("Vaccine-unique clusters found:")
print(vaccine_unique_df)

write.csv(vaccine_unique_df,
          file.path(OUTPUT_DIR, "vaccine_unique_clusters.csv"),
          row.names = FALSE)

# Plot: stacked bar of vaccine composition per cluster, highlight unique ones
vaccine_plot_list <- lapply(names(obj_list), function(lin) {
  md  <- get_meta(obj_list[[lin]], res_map[lin])
  
  tbl <- md %>%
    group_by(cluster, vaccine) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(pct = n_cells / sum(n_cells) * 100) %>%
    ungroup()
  
  unique_cls <- vaccine_unique_df %>%
    filter(lineage == lin) %>%
    pull(cluster)
  
  tbl$is_unique <- ifelse(tbl$cluster %in% unique_cls, "unique", "mixed")
  
  ggplot(tbl, aes(x = cluster, y = pct, fill = vaccine)) +
    geom_bar(stat = "identity") +
    geom_text(
      data = tbl %>% filter(is_unique == "unique") %>% distinct(cluster, is_unique),
      aes(x = cluster, y = 105, label = "*"),
      inherit.aes = FALSE, size = 5, color = "red"
    ) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = paste0(lin, " — Vaccine composition per cluster"),
      subtitle = "* = cluster 100% from one vaccine",
      x        = "Cluster",
      y        = "% of cells",
      fill     = "Vaccine"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

pdf(file.path(OUTPUT_DIR, "vaccine_unique_clusters.pdf"), width = 10, height = 5)
for (p in vaccine_plot_list) print(p)
dev.off()
message("Saved: vaccine_unique_clusters.pdf")

# =============================================================================
# Analysis 2 — Clusters appearing in only one patient
# =============================================================================

message("\n--- Analysis 2: Patient-specific clusters ---")

patient_unique_list <- lapply(names(obj_list), function(lin) {
  md  <- get_meta(obj_list[[lin]], res_map[lin])
  
  tbl <- md %>%
    group_by(cluster, patient_id) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(
      total_cells = sum(n_cells),
      n_patients  = n_distinct(patient_id),
      pct         = n_cells / total_cells * 100
    ) %>%
    ungroup()
  
  tbl %>%
    filter(n_patients == 1) %>%
    dplyr::select(cluster, patient_id, n_cells, total_cells, pct) %>%
    mutate(lineage = lin)
})

patient_unique_df <- bind_rows(patient_unique_list)

message("Patient-specific clusters found:")
print(patient_unique_df)

write.csv(patient_unique_df,
          file.path(OUTPUT_DIR, "patient_specific_clusters.csv"),
          row.names = FALSE)

# Plot: number of patients per cluster (bar), highlight single-patient ones
patient_plot_list <- lapply(names(obj_list), function(lin) {
  md  <- get_meta(obj_list[[lin]], res_map[lin])
  
  tbl <- md %>%
    group_by(cluster) %>%
    summarise(
      n_patients  = n_distinct(patient_id),
      total_cells = n(),
      .groups     = "drop"
    ) %>%
    mutate(is_unique = n_patients == 1)
  
  ggplot(tbl, aes(x = reorder(cluster, n_patients),
                  y = n_patients,
                  fill = is_unique)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "#6baed6", "TRUE" = "#fc4e2a"),
                      labels = c("FALSE" = "Multiple patients",
                                 "TRUE"  = "Single patient")) +
    labs(
      title = paste0(lin, " — Patients per cluster"),
      x     = "Cluster",
      y     = "N patients",
      fill  = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

pdf(file.path(OUTPUT_DIR, "patient_specific_clusters.pdf"), width = 10, height = 5)
for (p in patient_plot_list) print(p)
dev.off()
message("Saved: patient_specific_clusters.pdf")

# =============================================================================
# Analysis 3 — Cell count distribution per cluster within each lineage
# =============================================================================

message("\n--- Analysis 3: Cell count distribution per cluster ---")

count_dist_list <- lapply(names(obj_list), function(lin) {
  md  <- get_meta(obj_list[[lin]], res_map[lin])
  
  tbl <- md %>%
    group_by(cluster) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(
      pct     = n_cells / sum(n_cells) * 100,
      lineage = lin
    )
  tbl
})

count_dist_df <- bind_rows(count_dist_list)

write.csv(count_dist_df,
          file.path(OUTPUT_DIR, "cell_count_per_cluster.csv"),
          row.names = FALSE)

# Plot A: absolute counts per cluster, faceted by lineage
p_abs <- ggplot(count_dist_df,
                aes(x = reorder(cluster, -n_cells), y = n_cells, fill = lineage)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_hline(yintercept = 50,  linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey20", linewidth = 0.5) +
  facet_wrap(~lineage, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Cell counts per cluster",
    x     = "Cluster",
    y     = "N cells"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot B: percentage of lineage per cluster
p_pct <- ggplot(count_dist_df,
                aes(x = reorder(cluster, -pct), y = pct, fill = lineage)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~lineage, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Cell % per cluster (within lineage)",
    x     = "Cluster",
    y     = "% of lineage cells"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(OUTPUT_DIR, "cell_count_distribution.pdf"), width = 14, height = 8)
print(p_abs)
print(p_pct)
dev.off()
message("Saved: cell_count_distribution.pdf")


# =============================================================================
# Cluster filtering — remove noisy clusters
# Criteria (ALL must be met to keep):
#   - n_cells >= 50          (absolute floor)
#   - pct     >= 0.1         (>= 0.1% of lineage)
#   - not vaccine-unique     (100% cells from one vaccine)
#   - not patient-unique     (all cells from one patient)
# =============================================================================

message("\n--- Cluster filtering ---")

# Build a per-lineage per-cluster flag table
filter_flags <- count_dist_df %>%
  left_join(
    vaccine_unique_df %>%
      distinct(lineage, cluster) %>%
      mutate(is_vaccine_unique = TRUE),
    by = c("lineage", "cluster")
  ) %>%
  left_join(
    patient_unique_df %>%
      distinct(lineage, cluster) %>%
      mutate(is_patient_unique = TRUE),
    by = c("lineage", "cluster")
  ) %>%
  replace_na(list(is_vaccine_unique = FALSE, is_patient_unique = FALSE)) %>%
  mutate(
    keep_cluster = n_cells >= 50 & pct >= 0.1 & !is_vaccine_unique & !is_patient_unique
  )

message("Cluster filter summary:")
print(
  filter_flags %>%
    group_by(lineage) %>%
    summarise(
      total_clusters   = n(),
      clusters_kept    = sum(keep_cluster),
      clusters_removed = sum(!keep_cluster),
      .groups = "drop"
    )
)

message("\nRemoved clusters:")
print(
  filter_flags %>%
    filter(!keep_cluster) %>%
    dplyr::select(lineage, cluster, n_cells, pct, is_vaccine_unique, is_patient_unique) %>%
    arrange(lineage, cluster)
)

write.csv(filter_flags,
          file.path(OUTPUT_DIR, "cluster_filter_flags.csv"),
          row.names = FALSE)

# Apply filter: subset each Seurat object to kept clusters only
clusters_to_keep <- filter_flags %>%
  filter(keep_cluster) %>%
  dplyr::select(lineage, cluster)

obj_list_filtered <- lapply(names(obj_list), function(lin) {
  sobj     <- obj_list[[lin]]
  res_col  <- res_map[lin]
  keep_cls <- clusters_to_keep %>% filter(lineage == lin) %>% pull(cluster)
  cells_keep <- rownames(sobj@meta.data)[as.character(sobj@meta.data[[res_col]]) %in% keep_cls]
  subset(sobj, cells = cells_keep)
})
names(obj_list_filtered) <- names(obj_list)

# Report cell counts before/after
message("\nCell counts before/after filtering:")
before <- sapply(obj_list,          ncol)
after  <- sapply(obj_list_filtered, ncol)
count_comparison <- data.frame(
  lineage       = names(before),
  cells_before  = before,
  cells_after   = after,
  cells_removed = before - after,
  pct_removed   = round((before - after) / before * 100, 2)
)
print(count_comparison)
write.csv(count_comparison,
          file.path(OUTPUT_DIR, "cell_counts_before_after_filter.csv"),
          row.names = FALSE)

# Save filtered objects
saveRDS(obj_list_filtered,
        file.path(dirname(RDS_PATH), "subobjects_by_lineage_filtered.rds"))
message("Saved filtered objects: results/subobjects_by_lineage_filtered.rds")

# Plot: before/after cell counts per lineage
count_comparison_long <- count_comparison %>%
  select(lineage, cells_before, cells_after) %>%
  pivot_longer(cols = c(cells_before, cells_after),
               names_to  = "stage",
               values_to = "n_cells") %>%
  mutate(stage = factor(stage, levels = c("cells_before", "cells_after"),
                        labels = c("Before", "After")))

p_filter <- ggplot(count_comparison_long,
                   aes(x = lineage, y = n_cells, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Before" = "#74a9cf", "After" = "#0570b0")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Cell counts before and after cluster filtering",
    x     = "Lineage",
    y     = "N cells",
    fill  = NULL
  ) +
  theme_bw(base_size = 11)

pdf(file.path(OUTPUT_DIR, "filter_before_after.pdf"), width = 8, height = 5)
print(p_filter)
dev.off()
message("Saved: filter_before_after.pdf")

# =============================================================================
# Summary report (text)
# =============================================================================

sink(file.path(OUTPUT_DIR, "summary.txt"))
cat("=== Lineage Cluster Analysis Summary ===\n\n")

cat("Resolutions used:\n")
for (lin in names(res_map)) cat(sprintf("  %-10s %s\n", lin, res_map[lin]))

cat("\n--- Analysis 1: Vaccine-unique clusters (100% one vaccine) ---\n")
if (nrow(vaccine_unique_df) == 0) {
  cat("  None found.\n")
} else {
  print(vaccine_unique_df %>% dplyr::select(lineage, cluster, vaccine, n_cells, total_cells))
}

cat("\n--- Analysis 2: Patient-specific clusters (single patient_id) ---\n")
if (nrow(patient_unique_df) == 0) {
  cat("  None found.\n")
} else {
  print(patient_unique_df %>% dplyr::select(lineage, cluster, patient_id, n_cells, total_cells))
}

cat("\n--- Analysis 3: Cell counts per cluster ---\n")
print(count_dist_df %>% arrange(lineage, cluster))

cat("\n--- Cluster filtering summary ---\n")
print(count_comparison)
cat("\nFull filter flags saved to: cluster_filter_flags.csv\n")
cat("Filtered Seurat objects saved to: results/subobjects_by_lineage_filtered.rds\n")
sink()

message("\nDone. All outputs in: ", OUTPUT_DIR)
message("  vaccine_unique_clusters.csv / .pdf")
message("  patient_specific_clusters.csv / .pdf")
message("  cell_count_distribution.pdf")
message("  cell_count_per_cluster.csv")
message("  cluster_filter_flags.csv")
message("  cell_counts_before_after_filter.csv")
message("  filter_before_after.pdf")
message("  summary.txt")
message("  results/subobjects_by_lineage_filtered.rds")
