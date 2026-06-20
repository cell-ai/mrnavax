O output do cellranger foi baixado para cada estudo e foram utilizados thresholds para filtragem inicial das células (1_preprocessing_ind_studies)
```
seurat_qc <- subset(
  seurat_merged,
  subset = nCount_RNA       >= 500   &
           nCount_RNA       <= 60000 &
           nFeature_RNA     >= 500   &
           nFeature_RNA     <= 8000  &
           log10GenesPerUMI >= 0.80  &
           percent_mt       <  20
)
```
Os datasets foram integrados individualmente e depois feita a integração entre estudos (2_Integrate)

Fiz a clusterização e anotação inicial para achar grandes grupos celulares (res 0.5 - confirmar).
<img width="646" height="382" alt="image" src="https://github.com/user-attachments/assets/4be8fa8d-75e7-4b91-9bf1-df199fd674e9" />

Esses clusters foram utilizados para subsetar objectos, reclusterizar e morever clusters estranhos (4_Recluster_largetypes_n_reannotate)

Os barcodes das células mantidas foram utilizadas para filtrar o objeto principal e estou tentando novamente anotar os subclusters.
Para isso, eu fiz o AddModuleScore no objeto principal e separei baseado os clusters baseados no score para cada linhagem.

Fiz os embeddings com scgpt e estou reanotando.
