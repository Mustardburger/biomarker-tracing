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
options(future.globals.maxSize = 100000 * 1024^2)  # Increase to 1000GB

### A more efficient merge (http://thecodesearch.com/2023/03/19/solve-the-problem-of-excessive-memory-consumption-in-seurat-merge/)

### Manually change metadata in cortical brain atlases and tabula sapiens
base_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data"
brain = readRDS(paste0(base_path, "/human_cortical_cell_atlas/MTG_DFC_AnG_orig_comb_counts_downsamp_cell-type-20000.rds"))
sapiens = readRDS(paste0(base_path, "/tabula_sapiens/tabula_sapiens_all_cells_downsamp_cell-type-10000.rds"))

# All cols
cols = c("cell_type", "tissue", "CrossArea_subclass", "compartment", "free_annotation", "cell_ontology_class")

# Brain's cols: CrossArea_subclass, tissue, cell_type
sapiens_meta = sapiens[[]]
sapiens_meta$CrossArea_subclass = "NA"
sapiens_meta = sapiens_meta[, cols]
sapiens@meta.data = sapiens_meta

# Sapien's cols: tissue, cell_type, compartment, free_annotation, cell_ontology_class, 
brain_meta = brain[[]]
brain_meta$compartment = "brain"
brain_meta$free_annotation = "brain"
brain_meta$cell_ontology_class = "brain"
brain_meta = brain_meta[, cols]
brain@meta.data = brain_meta

# Merge the 2 metas
meta_total = rbind(brain_meta, sapiens_meta)


### Merge the two atlas
obj.ls <- list()
seu_l = c(sapiens, brain)
name_l = c("sapiens", "brain")

# Get the variable genes
var_genes = Features(seu_l[[1]])
for(i in 2:2){
    var = Features(seu_l[[i]])
    var_genes = union(var_genes, var)
}

# Extract dgcMatrix object from Seurat
# change the slot for data and counts for saving either normalized data or raw counts
for(i in 1:2){
    name = name_l[i]
    seu = seu_l[i][[1]]
    var = Features(seu)

    # colnames of dgCMatrix is cell barcodes
    obj.ls[[name]] <- GetAssayData(seu, layer = "counts", assay = "RNA")[var, ]
    print(paste0("Shape of ", name, ": ", dim(obj.ls[[name]])))
    #seu[["RNA"]]$counts
}

# Merge metadata by row, assuming all seurat has similar #columns
# If not, then have to process manually
metadata_comb = NULL
for(i in 1:2){
    seu = seu_l[i][[1]]
    if (is.null(metadata_comb)) {
        metadata_comb = seu[[]]
    } else {
        metadata_comb = rbind(metadata_comb, seu[[]])
    }
}

# Get union genes
get_all_genes <- function(obj.ls){
   all_genes <- rownames(obj.ls[[names(obj.ls)[1]]])
   for(f in names(obj.ls)[2:length(names(obj.ls))]){
      all_genes <- union(all_genes, rownames(obj.ls[[f]]))
   }
   return(all_genes)
}

# add 0
add_zero <- function(mt, all_genes){
   gene_left <- setdiff(all_genes, rownames(mt))
   left_mt <- as(matrix(0, ncol = ncol(mt), nrow = length(gene_left)), "dgCMatrix")
   colnames(left_mt)  <- colnames(mt)
   rownames(left_mt) <- gene_left
   mt <- rbind(mt, left_mt)
   return(mt[all_genes,,drop = F])
}

# Cbind matrix 
cbind_dgC_lst <- function(dgc_lst){
    merge_mt <- dgc_lst[[names(dgc_lst)[1]]]
    for(i in names(dgc_lst)[2:length(names(dgc_lst))]){
       merge_mt <- cbind(merge_mt, dgc_lst[[i]])
   }
   return(merge_mt)
}

# Get all genes
all_genes <- get_all_genes(obj.ls)
length(all_genes)

# Add zeroes to the matrix
new_obj.ls <- lapply(obj.ls, add_zero, all_genes)
rm(obj.ls)
gc()

# cbind the matrices
all_mt <- cbind_dgC_lst(new_obj.ls)
rm(new_obj.ls)
gc()
qsave(all_mt, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tmp/all_mt.qs")
write.table(metadata_comb, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tmp/metadata.tsv", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Save the final seurat obj
#all_mt = qread("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tmp/all_mt.qs")
#metadata_comb = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tmp/metadata.tsv", sep="\t")
combined <- CreateSeuratObject(counts = all_mt, min.cells = 100)
tmp = AddMetaData(combined, metadata = metadata_comb)
saveRDS(tmp, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens_cortical_cells/tabula_sapiens_cortical_cells_atlas_all-genes_downsamp.rds")

