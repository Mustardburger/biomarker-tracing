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

pseudobulk = AggregateExpression(obj, return.seurat=TRUE)
pseudo_exp_arr = pseudobulk[["RNA"]]$count
pseudo_exp_arr_log = pseudobulk[["RNA"]]$data
write.table(data.frame(pseudo_exp_arr), paste(save_path, "tabula_sapiens_pseudobulk_tissue_gene_exp_counts.tsv", sep="/"), row.names = TRUE, quote = FALSE, sep = "\t")
write.table(data.frame(pseudo_exp_arr_log), paste(save_path, "tabula_sapiens_pseudobulk_tissue_gene_exp_logcounts.tsv", sep="/"), row.names = TRUE, quote = FALSE, sep = "\t")

# Go to Jupyter notebook to plot out these gene markers
gene = "ENSG00000108849"
gene_row = sort(pseudo_exp_arr_log[which(rownames(pseudo_exp_arr_log) == gene), , drop = TRUE], decreasing = TRUE)
head(gene_row, 20)