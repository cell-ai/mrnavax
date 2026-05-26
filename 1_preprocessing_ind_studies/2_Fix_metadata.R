#!/usr/bin/env Rscript
# =============================================================
# 2_fix_metadata.R
# Extracts tissue / timepoint / donor from each Seurat object,
# compares them against the SRA run tables (ground truth),
# reports which samples are wrong, then saves corrected copies.
# Originals are NEVER overwritten.
#
# Output RDS:
#   results/GSE195673_seurat_corrected.rds
#   results/GSE328165_seurat_corrected.rds
#
# SRA columns used:
#   Sample Name  -> matches obj$sample
#   tissue       -> replaces obj$tissue   ("blood"->"PBMC", "lymph node"->"lymph_node")
#   timepoint    -> replaces obj$timepoint
#   isolate      -> replaces obj$donor
#   sex          -> added as new column
# =============================================================

setwd("/media/csbl/sandbox-SSD-3/mra_vaccine")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

dir.create("results/metadata_audit", recursive = TRUE, showWarnings = FALSE)

# =============================================================
# 1. Load SRA tables and build one-row-per-sample lookups
# =============================================================

# Multiple SRA runs exist per sample (GEX + BCR + TCR libraries)
# but tissue / timepoint / isolate are identical across runs,
# so we just deduplicate on Sample Name.

sra_195 <- read.csv("h5ad_geo/SraRunTable_GSE195673.csv",
                    stringsAsFactors = FALSE) %>%
  select(sample_id = Sample.Name, tissue, timepoint, donor = isolate, sex) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  mutate(
    tissue = case_when(
      tissue == "blood"       ~ "PBMC",
      tissue == "lymph node"  ~ "lymph_node",
      TRUE                    ~ tissue
    )
  )

sra_328 <- read.csv("h5ad_geo/SraRunTable_GSE328165.csv",
                    stringsAsFactors = FALSE) %>%
  select(sample_id = Sample.Name, tissue, timepoint, donor = isolate, sex) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  mutate(
    tissue = case_when(
      tissue == "blood"       ~ "PBMC",
      tissue == "lymph node"  ~ "lymph_node",
      TRUE                    ~ tissue
    )
  )

# =============================================================
# 2. Check and fix function
# =============================================================

fix_object <- function(study_id, rds_path, sra) {
  
  message("\n========== ", study_id, " ==========")
  obj  <- readRDS(rds_path)
  
  # ── Extract one row per sample from the object ───────────────
  # sex is new — pull it from @meta.data only if it already exists
  obj_meta <- obj@meta.data %>%
    select(sample, tissue, timepoint, donor,
           any_of("sex")) %>%
    distinct(sample, .keep_all = TRUE)
  
  check <- obj_meta %>%
    left_join(sra, by = c("sample" = "sample_id"),
              suffix = c("_obj", "_sra"))
  
  # Samples in object but missing from SRA table
  no_match <- filter(check, is.na(tissue_sra))
  if (nrow(no_match) > 0) {
    message("[WARN] ", nrow(no_match),
            " sample(s) not found in SRA table — left unchanged:")
    print(select(no_match, sample, tissue_obj, timepoint_obj, donor_obj),
          row.names = FALSE)
  }
  
  # ── Report mismatches / new columns ──────────────────────────
  wrong_tissue    <- filter(check, !is.na(tissue_sra),
                            tissue_obj    != tissue_sra)
  wrong_timepoint <- filter(check, !is.na(timepoint_sra),
                            timepoint_obj != timepoint_sra)
  wrong_donor     <- filter(check, !is.na(donor_sra),
                            donor_obj     != donor_sra)
  
  if (nrow(wrong_tissue) == 0 &&
      nrow(wrong_timepoint) == 0 &&
      nrow(wrong_donor) == 0) {
    message("[OK] tissue / timepoint / donor already match the SRA table.")
  } else {
    if (nrow(wrong_tissue) > 0) {
      message("\n[FIX] Tissue — ", nrow(wrong_tissue), " samples wrong:")
      print(select(wrong_tissue, sample, tissue_obj, tissue_sra),
            row.names = FALSE)
    }
    if (nrow(wrong_timepoint) > 0) {
      message("\n[FIX] Timepoint — ", nrow(wrong_timepoint), " samples wrong:")
      print(select(wrong_timepoint, sample, timepoint_obj, timepoint_sra),
            row.names = FALSE)
    }
    if (nrow(wrong_donor) > 0) {
      message("\n[FIX] Donor — ", nrow(wrong_donor), " samples wrong:")
      print(select(wrong_donor, sample, donor_obj, donor_sra),
            row.names = FALSE)
    }
  }
  
  message("\n[NEW] Adding sex column from SRA table.")
  
  # ── Apply corrections + add sex directly to @meta.data ───────
  lut <- sra %>% filter(sample_id %in% obj@meta.data$sample)
  
  lut_tissue    <- setNames(lut$tissue,    lut$sample_id)
  lut_timepoint <- setNames(lut$timepoint, lut$sample_id)
  lut_donor     <- setNames(lut$donor,     lut$sample_id)
  lut_sex       <- setNames(lut$sex,       lut$sample_id)
  
  s         <- obj@meta.data$sample
  has_match <- s %in% lut$sample_id
  
  obj@meta.data$tissue[has_match]    <- lut_tissue   [s[has_match]]
  obj@meta.data$timepoint[has_match] <- lut_timepoint[s[has_match]]
  obj@meta.data$donor[has_match]     <- lut_donor    [s[has_match]]
  obj@meta.data$sex                  <- NA_character_
  obj@meta.data$sex[has_match]       <- lut_sex      [s[has_match]]
  
  # ── Verify ────────────────────────────────────────────────────
  message("\nPost-fix check:")
  for (col in c("tissue", "timepoint", "donor", "sex")) {
    vals <- unique(obj@meta.data[[col]])
    message("  ", col, ": ", paste(sort(vals), collapse = ", "))
  }
  
  # ── Save (new file — original is never touched) ───────────────
  out_rds <- sub("_seurat_processed\\.rds$", "_seurat_corrected.rds", rds_path)
  message("\nSaving corrected object: ", out_rds)
  message("Original preserved at  : ", rds_path)
  saveRDS(obj, out_rds)
  
  out_csv <- paste0("results/metadata_audit/", study_id, "_metadata_fixed.csv")
  obj@meta.data %>%
    group_by(sample, gse, tissue, timepoint, donor, sex) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    arrange(donor, timepoint) %>%
    write.csv(out_csv, row.names = FALSE)
  message("Summary saved: ", out_csv)
  
  invisible(obj)
}

# =============================================================
# 3. Run
# =============================================================


fix_object("GSE328165",
           "results/GSE328165_seurat_processed.rds",
           sra_328)

fix_object("GSE195673",
           "results/GSE195673_seurat_processed.rds",
           sra_195)

message("\nDone.")
