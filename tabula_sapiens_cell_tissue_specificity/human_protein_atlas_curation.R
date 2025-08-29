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
  library("glue")
  library(stringr)
  library(purrr)
})
options(future.globals.maxSize = 10000 * 1024^2)


# Read in the main file
cluster_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/human_protein_atlas_single_cell"
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/human_protein_atlas_single_cell/cell_tissue/specificity_metric"
clus_df = fread(glue("{cluster_path}/rna_single_cell_cluster.tsv"))

# Cell tissue
pseudobulk_cell = clus_df %>% 
    rename(cell_type = "Cell type", read_count = "Read count", tissue = "Tissue", gene = "Gene") %>% 
    group_by(tissue, cell_type, gene) %>% summarise(total_read_count = sum(read_count), .groups = "drop") %>%
    mutate(cell_type = gsub(" |-", ".", cell_type), tissue = gsub(" |-", ".", tissue)) %>%
    mutate(cell_tissue = paste(tissue, cell_type, sep = "_")) %>%
    select(gene, cell_tissue, total_read_count) %>%
    pivot_wider(names_from = cell_tissue, values_from = total_read_count, values_fill = 0)

# Loop over tissues and see if there is any cluster with fewer than 50 cells
# Use 2 nested for loops here, should be better solution
files = Sys.glob(glue("{cluster_path}/*/cell_data.tsv"))
cells_to_drop = c(); thres = 50
for (file in files) {
    tissue = strsplit(file, split = "/")[[1]][13]
    print(glue("Working on {tissue}"))
    a = fread(glue("{cluster_path}/{tissue}/cell_data.tsv"))

    clus_to_cell_type = clus_df %>% 
        dplyr::mutate(clus_id = map_chr(str_split(Cluster, "-"), 2)) %>%
        dplyr::filter(Tissue == gsub("_", " ", tissue)) %>% 
        dplyr::select(`Cell type`, clus_id) %>%
        dplyr::group_by(`Cell type`) %>%
        dplyr::summarise(uniq_clus = list(unique(clus_id)))
    if (nrow(clus_to_cell_type) == 0){
      break
    }
    
    for (i in 1:nrow(clus_to_cell_type)) {
      cell_type = gsub(" |-", ".", clus_to_cell_type[["Cell type"]][i]); tissue_name = gsub(" |_", ".", tissue)
      print(cell_type)
      cell_tissue = paste(tissue_name, cell_type, sep = "_")
      l = clus_to_cell_type[["uniq_clus"]][i] %>% unlist() %>% as.numeric()
      cell_count = a[a$cluster %in% l, ] %>% nrow()

      if (cell_count < thres) {
        cells_to_drop = c(cells_to_drop, cell_tissue)
      }
    }
}
cells_to_drop = c("bone.marrow_plasma.cells", "esophagus_schwann.cells", "pancreas_endothelial.cells", "pancreas_macrophages", "pbmc_dendritic.cells", "small.intestine_mixed.immune.cells", "testis_sertoli.cells")
pseudobulk_cell = pseudobulk_cell %>% as.data.frame() %>% select(-cells_to_drop)

fwrite(pseudobulk_cell, glue("{save_path}/human_protein_atlas_sc_pseudobulk_gene_exp_counts.tsv"),
  sep = "\t", col.names = T
)

# Convert to a Seurat object
rownames(pseudobulk_cell) <- pseudobulk_cell$gene
pseudobulk_cell$gene <- NULL
expr_matrix <- as.matrix(pseudobulk_cell); rownames(expr_matrix) = rownames(pseudobulk_cell)
seurat_obj <- CreateSeuratObject(counts = expr_matrix)
seurat_obj = NormalizeData(seurat_obj)
pseudo_exp_arr_log = seurat_obj[["RNA"]]$data %>% as.data.frame()
pseudo_exp_arr_log$gene = rownames(pseudo_exp_arr_log)
rownames(pseudo_exp_arr_log) = NULL

fwrite(pseudo_exp_arr_log, glue("{save_path}/human_protein_atlas_sc_pseudobulk_gene_exp_logcounts.tsv"),
  sep = "\t", col.names = T
)
saveRDS(seurat_obj, glue("{save_path}/human_protein_atlas_sc_pseudobulk_seurat.rds"))


### Tissue level
cell_tissue = rownames(seurat_obj[[]])
tissue = lapply(cell_tissue, function (x) strsplit(x, split = "_")[[1]][1]) %>% unlist()
seurat_obj[[]]$tissue = tissue
Idents(seurat_obj) = "tissue"
pseudobulk_tis = AggregateExpression(seurat_obj, return.seurat=TRUE)
pseudo_exp_arr_tis = pseudobulk_tis[["RNA"]]$count %>% as.data.frame()
pseudo_exp_arr_log_tis = pseudobulk_tis[["RNA"]]$data %>% as.data.frame()

pseudo_exp_arr_tis$gene = rownames(pseudo_exp_arr_tis)
rownames(pseudo_exp_arr_tis) = NULL
pseudo_exp_arr_log_tis$gene = rownames(pseudo_exp_arr_log_tis)
rownames(pseudo_exp_arr_log_tis) = NULL
fwrite(pseudo_exp_arr_tis, glue("{save_path}/human_protein_atlas_sc_pseudobulk_tissue_gene_exp_counts.tsv"),
  sep = "\t", col.names = T
)
fwrite(pseudo_exp_arr_log_tis, glue("{save_path}/human_protein_atlas_sc_pseudobulk_tissue_gene_exp_logcounts.tsv"),
  sep = "\t", col.names = T
)

# Some sanity check
gene_name = "ENSG00000277586"
pseudo_exp_arr_log %>% 
    dplyr::filter(gene == gene_name) %>%
    pivot_longer(cols = -gene) %>% arrange(desc(value))