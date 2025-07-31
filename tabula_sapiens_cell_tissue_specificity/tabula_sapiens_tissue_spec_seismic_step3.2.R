suppressMessages({
  library("plyr")
  library("ggrepel")
  library("ggthemes")
  library("grid")
  library("Seurat")
  #library("SeuratDisk")
  library("SeuratData")
  library("scales")
  library("qs")
  library("scRNAseq")
  library("dplyr")
  library("tidyr")
  library("data.table")
  library("ggplot2")
})
options(future.globals.maxSize = 10000 * 1024^2)


### Step 3.1: Calculate cell type specificity based on average gene expression of tissue-cell-type pair
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/tabula_sapiens/cell_tissue/specificity_metric"
dir.create(save_path, showWarnings = FALSE)

# Read in the data (461 unique cell-type-tissue pairs)
obj = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_50_above_cells_per_cell_tissue.rds")
Idents(obj) = "tissue"


library("seismicGWAS")
DefaultAssay(obj) = "RNA"

obj_diet = DietSeurat(obj, assays = "RNA", layers = "data")
obj_sce = as.SingleCellExperiment(obj_diet)
obj_sscore <- calc_specificity(obj_sce, ct_label_col='tissue', 
                             min_avg_exp_ct=0.0, min_uniq_ct = 0,
                             min_ct_size = 0, min_cells_gene_exp = 0
                             )

seis_df = as.data.frame(obj_sscore)
colnames(seis_df) <- sapply(colnames(seis_df), function(i) {
  parts <- strsplit(i, ":")[[1]]
  paste0(gsub(" ", ".", parts[2]), "_", gsub(" ", ".", parts[1]))
})
file_path = glue("{save_path}/tabula_sapiens_tissue_seismic.tsv")
write.table(as.data.frame(seis_df), file = file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Filter only genes in Olink
atlas_data_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/atlas_data"
atlas_all_cells = fread(glue("{atlas_data_path}/atlas_all_cell_tissues.tsv"))
seis_df_smal = seis_df[olink_genes, ]
seis_df_smal$gene = rownames(seis_df_smal)
rownames(seis_df_smal) = NULL
fwrite(seis_df_smal, glue("{atlas_data_path}/atlas_seismic_all_cell_tissues.tsv"), row.names = F, col.names = T, sep = "\t")