setwd("/mnt/sandbox-SSD/marcela_ishihara/project/mravax/")

library(DESeq2)
library(volcano3D)
library(Seurat)
library(plotly)
library(DESeq2)

####adding metadata

seurat_obj <- readRDS("merged_after_reclustering_w_res_dims20.rds")
meta <- seurat_obj[[]]

sample_meta <- meta |>
  dplyr::distinct(sample, sex, timepoint, study)

View(sample_meta)

library(dplyr)
library(ggplot2)

meta <- seurat_obj[[]]


library(dplyr)
library(Seurat)

sra1 <- read.csv("SraRunTable (4).csv", check.names = FALSE)
sra2 <- read.csv("SraRunTable (5).csv", check.names = FALSE)

sra <- bind_rows(sra1, sra2)

library_df <- sra %>%
  filter(`Assay Type` == "RNA-Seq") %>%
  transmute(
    sample = `Sample Name`,
    library_name = `Library Name`,
    library_run = Run,
    library_biorep = Biological_replicate,
    sra_study = `SRA Study`
  ) %>%
  distinct()


library_df %>%
  count(sample) %>%
  filter(n > 1)

meta <- seurat_obj[[]]

meta2 <- meta %>%
  tibble::rownames_to_column("cell") %>%
  left_join(library_df, by = "sample") %>%
  tibble::column_to_rownames("cell")

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = meta2[, c("library_name", "library_run", "library_biorep", "sra_study")]
)

head(seurat_obj@meta.data[, c("sample", "library_name", "library_run")])

library(dplyr)
library(stringr)
library(Seurat)

meta <- seurat_obj[[]] %>%
  mutate(
    patient_id = stringr::str_extract(library_name, "WU397-\\d{3}|WU368-\\d{2}"),
    vaccine = case_when(
      stringr::str_starts(library_name, "ELAB-WU368") ~ "covidmRNA",
      patient_id %in% c("WU397-006", "WU397-017", "WU397-028", "WU397-029") ~ "mRNA-1010",
      patient_id %in% c("WU397-005", "WU397-009", "WU397-022") ~ "Fluarix",
      TRUE ~ NA_character_
    )
  )

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = meta[, c("patient_id", "vaccine")]
)

seurat_obj[[]] %>%
  distinct(library_name, patient_id, vaccine) %>%
  arrange(patient_id) %>%
  View()

table(seurat_obj$vaccine, useNA = "ifany")

### compositional plots
library(dplyr)
library(ggplot2)

meta <- seurat_obj[[]]

plot_cluster_composition <- function(meta, variable, cluster_col = "RNA_snn_res.0.1") {
  
  cluster_meta <- meta %>%
    count(
      cluster = .data[[cluster_col]],
      variable = .data[[variable]],
      name = "n_cells"
    )
  
  ggplot(cluster_meta, aes(x = cluster, y = n_cells, fill = variable)) +
    geom_col(position = "fill") +
    theme_bw() +
    labs(
      x = "Cluster",
      y = "Proportion of cells",
      fill = variable,
      title = paste("Cluster composition by", variable)
    )
}

plot_cluster_composition(meta, "study")
plot_cluster_composition(meta, "sex")
plot_cluster_composition(meta, "vaccine")
plot_cluster_composition(meta, "patient_id")
plot_cluster_composition(meta, "timepoint")



######three way volcano
counts_mat <-AggregateExpression(seurat_obj, 
                                 group.by = c("vaccine", "RNA_snn_res.0.1",
                                              "timepoint"),
                                 assays = "RNA",
                                 slot = "counts",
                                 return.seurat = FALSE)

raw_counts <-as.data.frame(counts_mat$RNA)

# Extract tissue and timepoint from the column names
vaccine <- gsub("_.*", "", colnames(raw_counts))  # Everything before the first underscore
timepoint <- gsub(".*_", "", colnames(raw_counts))  # Everything after the last underscore

# Create a data frame with these columns
metadata_test <- data.frame(colnames(raw_counts),
                            vaccine = vaccine, timepoint = timepoint)

colnames(metadata_test) <- c("sample", "vaccine", "timepoint")

# Make sure the tissue column exists and is a factor with 3 levels
metadata_test$vaccine <- factor(metadata_test$vaccine,
                               levels = c("Fluarix", "mRNA-1010", "covidmRNA"))

# Create DESeqDataSet using tissue as the design variable
dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = metadata_test, 
                              design = ~ vaccine)

dds_DE <- DESeq(dds)

dds_LRT <- DESeq(dds, test = "LRT", reduced = ~ 1, parallel = TRUE)

res <- deseq_polar(dds_DE,dds_LRT, "vaccine", padj.method = "fdr")

radial_ggplot(polar = res,
              marker_size = 2.3,
              marker_outline_width = 0,
              legend_size = 10,
              plot_top = 40) +        # <-- limits the -log10(p) axis to 40
  theme(legend.position = "right")

# Plot the 3D volcano plot
volcano3D(res)

# -log10(x) = 40  means  x = 10^-40
# So cap padj at a minimum of 10^-40 (prevents z going above 40)
# Check the df slot — it's a LIST of 2 dataframes
str(slot(res, "df"))

# Look at the z column in each dataframe
head(slot(res, "df")[[1]])   # type 1 (scaled)
head(slot(res, "df")[[2]])   # type 2 (unscaled)
volcano3D(res)



########################################################
early <- c("d0", "d8", "d15")
mid   <- c("d28", "d28+d35", "d35", "d57", "d60")
late  <- c("d110", "d121", "d180", "d181", "d201")

metadata_test$timegroup <- case_when(
  metadata_test$timepoint %in% early ~ "early",
  metadata_test$timepoint %in% mid   ~ "mid",
  metadata_test$timepoint %in% late  ~ "late"
)

metadata_test$timegroup <- factor(metadata_test$timegroup,
                                  levels = c("early", "mid", "late"))

# Verify the new distribution
table(metadata_test$vaccine, metadata_test$timegroup)


###############################helder
counts_mat <- AggregateExpression(seurat_obj, 
                                  group.by = c("vaccine", "RNA_snn_res.0.1",
                                               "timepoint", "sex"),
                                  assays = "RNA",
                                  slot = "counts",
                                  return.seurat = FALSE)

raw_counts <- as.data.frame(counts_mat$RNA)

# Rename covidmRNA to BNT162b2 in column names
colnames(raw_counts) <- gsub("covidmRNA", "BNT162b2", colnames(raw_counts))

# Extract vaccine and timepoint from column names
vaccine   <- gsub("_.*", "", colnames(raw_counts))
timepoint <- gsub(".*_", "", colnames(raw_counts))

metadata_test <- data.frame(colnames(raw_counts),
                            vaccine   = vaccine,
                            timepoint = timepoint)
colnames(metadata_test) <- c("sample", "vaccine", "timepoint")

metadata_test$vaccine <- factor(metadata_test$vaccine,
                                levels = c("Fluarix", "mRNA-1010", "BNT162b2"))

dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData   = metadata_test, 
                              design    = ~ vaccine)

dds_DE  <- DESeq(dds)
dds_LRT <- DESeq(dds, test = "LRT", reduced = ~ 1, parallel = TRUE)

res <- deseq_polar(dds_DE, dds_LRT, "vaccine", padj.method = "fdr")

# Cap z axis at 40
slot(res, "df")[[1]]$z <- pmin(slot(res, "df")[[1]]$z, 40)
slot(res, "df")[[2]]$z <- pmin(slot(res, "df")[[2]]$z, 40)

radial_ggplot(polar                = res,
              marker_size          = 2.3,
              marker_outline_width = 0,
              legend_size          = 10,
              plot_top             = 40) +
  theme(legend.position = "right")

volcano3D(res)
library(htmlwidgets)

p <- volcano3D(res)

saveWidget(p, file = "volcano3D.html", selfcontained = TRUE)

