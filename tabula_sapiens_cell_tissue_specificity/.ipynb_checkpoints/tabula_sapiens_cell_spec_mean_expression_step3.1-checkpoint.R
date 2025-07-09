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
Idents(obj) = "cell_type_tissue"
obj[[]] = obj[[]] %>% mutate(cell_type_tissue2 = gsub("_", ":", cell_type_tissue))
Idents(obj) = "cell_type_tissue2"
#length(unique(obj[[]]$cell_type_tissue))

# Calculate the Average feature expression (strangely if layer=="data", the function automatically exponential the matrix first)
# avg_exp = AverageExpression(obj, layer="count", return.seurat=TRUE)
# avg_exp_arr = avg_exp[["RNA"]]$data
# cell_tissue = colnames(avg_exp_arr)

# Calculate the Pseudobulk expression
pseudobulk = AggregateExpression(obj, return.seurat=TRUE)
pseudo_exp_arr = pseudobulk[["RNA"]]$count
pseudo_exp_arr_log = pseudobulk[["RNA"]]$data
write.table(data.frame(pseudo_exp_arr), paste(save_path, "tabula_sapiens_pseudobulk_gene_exp_counts.tsv", sep="/"), row.names = TRUE, quote = FALSE, sep = "\t")
write.table(data.frame(pseudo_exp_arr_log), paste(save_path, "tabula_sapiens_pseudobulk_gene_exp_logcounts.tsv", sep="/"), row.names = TRUE, quote = FALSE, sep = "\t")

# Sanity check: look at individual genes (kind of correct)
gene = "ENSG00000185915"
gene_row = sort(avg_exp_arr[which(rownames(avg_exp_arr) == gene), , drop = TRUE], decreasing = TRUE)
head(gene_row)

# Sanity check: look at individual cell-tissue (looks correct)
cel_tis = "liver_hepatocyte"
cell_tissue_df = data.frame(cell_tissue = cell_tissue) %>% separate(cell_tissue, into=c("cell", "tissue"), sep="_")
gene_col = sort(avg_exp_arr[ , which(colnames(avg_exp_arr) == cel_tis), drop = TRUE], decreasing = TRUE)

# Save the expression to the df for easier processing
save_df = data.frame(avg_exp_arr)
write.table(save_df, paste(save_path, "tabula_sapiens_log_avg_gene_exp.tsv", sep="/"), row.names = TRUE, quote = FALSE, sep = "\t")


### Addendum: Calculate variable features and match them with Olink proteins
# Load in stuff
N = 29000
pseudo_var_obj = FindVariableFeatures(pseudobulk, nfeatures=N)
var_features = VariableFeatures(pseudo_var_obj)


## Save variable genes with union with Olink proteins
plasma_proteome_gene_df = read.table(
  file = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data/plasma_proteome_gene_id_name_map.tsv",
  sep = "\t", header = TRUE
)

# Intersect the two sets
common_vals = intersect(var_features, plasma_proteome_gene_df$gene)
length(common_vals)
all_vals = union(var_features, plasma_proteome_gene_df$gene)
length(all_vals)
labels = all_vals %in% var_features

# Save results
df = data.frame(gene = all_vals, label = labels)
write.table(df, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data/var_genes.tsv", quote=FALSE, row.names=FALSE, sep="\t")



## Save variable genes irrespective of Olink proteins
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/tabula_sapiens/cell_tissue/variable_features"
png(paste0(save_path, "/variable_feature_plot.png"), width=800, height=800)
VariableFeaturePlot(pseudo_var_obj)
dev.off()

# Save the matrix that contains ranks, mean, and standardized variance for each of the features
var_features_metadata = pseudo_var_obj[["RNA"]]@meta.data
write.table(var_features_metadata, paste0(save_path, "/var_genes.tsv"), quote=FALSE, row.names=FALSE, sep="\t")
