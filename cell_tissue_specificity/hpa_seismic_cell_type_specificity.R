# Install these packages first before running the code
library("seismicGWAS")
library("plyr")
library("Seurat")
library("qs")
library("scRNAseq")
library("dplyr")
library("data.table")
options(future.globals.maxSize = 10000 * 1024^2)


######################### RUNNING SEISMIC CELL TYPE SPECIFICITY #########################
# Read in the obj_smal
obj_smal = readRDS("tutorial/human_protein_atlas_single_cell_small.rds")

# Run cell type specificity
DefaultAssay(obj_smal) = "RNA"
obj_diet = DietSeurat(obj_smal, assays = "RNA", layers = "data")
obj_sce = as.SingleCellExperiment(obj_diet)
obj_sscore <- calc_specificity(
    obj_sce, ct_label_col=ident,
    min_uniq_ct = 0, min_ct_size = 0, min_cells_gene_exp = 0, min_avg_exp_ct = 0.0
)

# Save the result
df_ss = as.data.frame(obj_sscore)
df_ss$gene = rownames(df_ss)
rownames(df_ss) = NULL
write.table(
    df_ss, 
    file = "tutorial/human_protein_atlas_seismic_cell_spec_matrix.tsv", 
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)