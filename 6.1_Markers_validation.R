###the obj is merged_after_reclustering_no_TCR_BCR.rds

####cluster 0 - naive / central memory CD4 T
cd4_markers <- c("CD27", "CCR7", "SELL", "TCF7", "LEF1", "ITGB1") #ITGB1 should not appear

FeaturePlot(
  seurat_obj,
  features = cd4_markers,
  ncol = 4
)

DotPlot(
  seurat_obj,
  features = cd4_markers,
  group.by = "RNA_snn_res.0.1"
) +
  RotatedAxis()

####cluster 1 - activated/memory b cell
memb_markers <- c("CD27", "AIM2", "MS4A1", "FCRL4", "FCRL5")
FeaturePlot(
  seurat_obj,
  features = memb_markers,
  ncol = 4
)

DotPlot(
  seurat_obj,
  features = memb_markers,
  group.by = "RNA_snn_res.0.1"
) +
  RotatedAxis()

####cluster 2 - naive b cell | IL4R and FCER2 genes similar to Transitional B cells, but have lower expression of MME, CD24, and CD9
naiveb_markers <- c("IL4R", "FCER2", "MME", "CD24", "CD9")
FeaturePlot(
  seurat_obj,
  features = naiveb_markers,
  ncol = 4
)

####cluster 3 - plasma b cells
plasma_markers <- c("CD27", "CD38", "PRDM1", "XBP1", "MZB1", "SLAMF7")
FeaturePlot(
  seurat_obj,
  features = plasma_markers,
  ncol = 4
)

####cluster 4 - cytotoxic cd8 t cells
cd8_markers <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B")
FeaturePlot(
  seurat_obj,
  features = cd8_markers,
  ncol = 5
)

####cluster 5 - GC B cells
#gcb_markers <- c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4")
gc_markers <- c("BCL6", "AICDA", "MEF2B", "S1PR2", "RGS13", "CD83", "CXCR4", "MKI67", "HMMR", "EZH2") #chatgpt
FeaturePlot(
  seurat_obj,
  features = gc_markers,
  ncol = 4
)

####cluster 6 - nk cells
nk_markers <- c("NKG7", "KLRC1", "NCAM1", "FCER1G", "FCGR3A", "GZMB") #allen e chatgpt
FeaturePlot(
  seurat_obj,
  features = nk_markers,
  ncol = 4
)

####cluster 7 - tfh cells
tfh_markers <- c("CXCR5", "PDCD1", "ICOS", "IL21", "MAF", "SH2D1A", "CD40LG") #chatgpt
FeaturePlot(
  seurat_obj,
  features = tfh_markers,
  ncol = 5
)

####cluster 8 - classical monocytes
mon_markers <- c("S100A8", "S100A9", "NFKBIA", "CD14", "VCAN") 
FeaturePlot(
  seurat_obj,
  features = mon_markers,
  ncol = 5
)

####cluster 9 - activated/memory T
memt_markers <- c("IL7R", "CD69", "ICOS", "LTB", "TNFRSF4", "TNFRSF9") #allen e chatgpt
FeaturePlot(
  seurat_obj,
  features = memt_markers,
  ncol = 5
)

####cluster 10 - plasmablasts
blast_markers <- c("CD27", "SLAMF7", "PRDM1", "MKI67", "TOP2A", "STMN1", "HMGB2", "JCHAIN") #chatgpt
FeaturePlot(
  seurat_obj,
  features = blast_markers,
  ncol = 5
)

####cluster 11 - pDC
pdc_markers <- c("PLAC8", "IRF8", "IL3RA", "ITM2C")
FeaturePlot(
  seurat_obj,
  features = pdc_markers,
  ncol = 5
)
