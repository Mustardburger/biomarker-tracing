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


# Open the file
seu = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells.rds")
seu[[]]$tissue_cell = paste0(seu[[]]$tissue, ":", seu[[]]$cell_type)
table(seu[[]]$tissue)
table(seu[[]]$cell_type)
median(table(seu[[]]$tissue_cell))

### Subset dataset so that each da
N = 500
full_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells.rds"
dataset_name = strsplit(basename(full_path), split="[.]")[[1]][1]

# Downsample code from Claude
ident = "tissue_cell"
Idents(seu) <- ident
cell_types <- unique(seu[[]][ident])[[1]]
downsampled_cells <- list()

for (ct in cell_types) {
  # Get cell barcodes of current cell type
  cells_of_type <- WhichCells(seu, idents = ct)
  
  # Sample cells from each cell type
  if (length(cells_of_type) > N) {
    sampled_cells <- sample(cells_of_type, N)
  } else {
    sampled_cells <- cells_of_type
  }
  
  downsampled_cells[[ct]] <- sampled_cells
}

# Subset object
all_downsampled_cells <- unlist(downsampled_cells)
seu_sub <- subset(seu, cells = all_downsampled_cells)

# Print out stuff
print(paste0(">>> ", dataset_name, " orig Seurat object has ", dim(seu)[2], " cells"))
print(paste0(">>> ", dataset_name, " subset Seurat object has ", dim(seu_sub)[2], " cells"))

# Save subsetted dataset
saveRDS(seu_sub, paste0(dirname(full_path), "/", dataset_name, "_subset.rds"))


### Do some basic quality
seu_orig = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells.rds")
seu_sub = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_subset.rds")

seu_orig[["percent.mt"]] = PercentageFeatureSet(seu_orig, pattern = "^MT-")


### Count number of cells
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/quality_control/tabula_sapiens_all_cells"
seu_type = "orig"
seu_orig[[]]$tissue_cell = paste0(seu_orig[[]]$tissue, ":", seu_orig[[]]$cell_type)
counts_c = table(seu_orig[[]]$tissue_cell)
#seu_sub[[]]$tissue_cell = paste0(seu_sub[[]]$tissue, ":", seu_sub[[]]$cell_type)
#counts_c = table(seu_sub[[]]$tissue_cell)

count_df = data.frame(counts_c, names(counts_c))

# Make count df
count_df = count_df %>% 
    select(!Var1) %>% 
    rename("count" = "Freq") %>%
    rename("tissue.cell" = "names.counts_c.") %>%
    separate(tissue.cell, into = c("tissue", "cell"), sep = ":")

# Make boxplot for tissue
p = ggplot(count_df, aes(x=reorder(tissue, count, FUN=median), y=count)) +
    scale_y_continuous(trans='log10') +
    xlab("tissue") +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
png(paste0(save_path, "/tissue_cell_count_", seu_type, ".png"), width = 1000, height = 800)
print(p)
dev.off()

# Make boxplot for cell type
p = ggplot(count_df, aes(x=reorder(cell, count, FUN=median), y=count)) +
    scale_y_continuous(trans='log10') +
    xlab("cell") +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
png(paste0(save_path, "/cell_type_cell_count_", seu_type, ".png"), width = 1200, height = 600, res=100)
print(p)
dev.off()

# Distribution of total cell counts
p1 = ggplot(count_df, aes(x=count)) + 
    geom_histogram()
p2 = ggplot(count_df, aes(y=count)) + 
    geom_boxplot()
p = (p1 | p2) + ggtitle(paste0("Cell ", seu_type))
png(paste0(save_path, "/cell_count_dist_", seu_type, ".png"), width = 1200, height = 600, res=300)
print(p)
dev.off()


### Remove genes with fewer than 3 cells
seu_orig = seu
mat_counts = seu_orig@assays$RNA@counts
mat_data =  seu_orig@assays$RNA@data
metadata = seu_orig[[]]

count_list = list(mat_counts, mat_data)
names(count_list) = c("counts", "data")
obj = CreateSeuratObject(counts = mat_counts, min.cells = 100, meta.data = metadata)
saveRDS(obj, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_fewer_genes.rds")


### Remove cell types with fewer than 100 cells
obj = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_fewer_genes.rds")
min_cells = 100
obj[[]]$tissue_cell = paste0(obj[[]]$tissue, ":", obj[[]]$cell_type)
tissue_cell_count = table(obj[[]]$tissue_cell)

tissue_cell_count_above_min = Filter(function(x) x > min_cells, tissue_cell_count)
tissue_cell_relevant = names(tissue_cell_count_above_min)

a = subset(x = obj, subset = tissue_cell %in% tissue_cell_relevant)
saveRDS(a, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_fewer_genes.rds")


### Retain only top variable genes
num_hvgs = 25000
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = num_hvgs)
saveRDS(obj, paste0("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_hvg", num_hvgs, ".rds"))


### Downsample number of cells for each cell type (regardless of tissue)
obj = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_fewer_genes.rds")
N = 10000
Idents(obj) = "cell_type"
obj_sub = subset(x = obj, downsample = N)
saveRDS(obj_sub, paste0("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_downsamp_cell-type-", N, ".rds"))


### Debug
tmp = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_fewer_genes.rds")
tmp2 = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells_hvg.rds")
g = Features(tmp)
g2 = VariableFeatures(obj)
gene_exa = "ENSG00000122711"
gene_exa = "ENSG00000099834"
gene_exa %in% g2