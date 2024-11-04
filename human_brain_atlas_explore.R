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
  library("scCustomize")
  library("ggplot2")
})

# CELLxGENE says the Seurat objects are v4
# Let's try if Seurat works here
base_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data"
#seu_lec = readRDS(paste0(base_path, "/human_brain_cell_atlas/Cx-AG-LEC.rds"))
#seu_mec = readRDS(paste0(base_path, "/human_brain_cell_atlas/Cx-APH-MEC.rds"))
#seu_mtg = readRDS(paste0(base_path, "/human_cortical_cell_atlas/MTG.rds"))
#seu_dfc = readRDS(paste0(base_path, "/human_cortical_cell_atlas/DFC.rds"))
#seu_ang = readRDS(paste0(base_path, "/human_cortical_cell_atlas/AnG.rds"))

N = 1000
paths = c(
  "/human_brain_cell_atlas/Cx-AG-LEC.rds",
  "/human_brain_cell_atlas/Cx-APH-MEC.rds",
)

# For brain_cell_atlas
for (path in paths) {
  full_path = paste0(base_path, path)
  dataset_name = strsplit(basename(full_path), split="[.]")[[1]][1]
  print(paste0("Working on ", dataset_name))
  seu = readRDS(full_path)

  # Downsample code from Claude
  ident = "supercluster_term"
  Idents(seu) <- ident
  cell_types <- unique(seu[[]][ident])
  downsampled_cells <- list()

  for (ct in cell_types) {
    # Get cell barcodes of current cell type
    print(ct)
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

}


paths = c(
  "/human_cortical_cell_atlas/MTG.rds",
  "/human_cortical_cell_atlas/DFC.rds",
  "/human_cortical_cell_atlas/AnG.rds"
)
# For cortical cell atlas
for (path in paths) {
  full_path = paste0(base_path, path)
  dataset_name = strsplit(basename(full_path), split="[.]")[[1]][1]
  print(paste0("Working on ", dataset_name))
  seu = readRDS(full_path)

  # Downsample code from Claude
  ident = "cell_type"
  Idents(seu) <- ident
  cell_types <- unique(seu[[]][ident])$cell_type
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

}

### A more efficient merge (http://thecodesearch.com/2023/03/19/solve-the-problem-of-excessive-memory-consumption-in-seurat-merge/)
# But in the end we still end up with a huge scRNA-seq dataset

### Count cells
seu_orig = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/human_cortical_cell_atlas/MTG_DFC_AnG_orig_comb_data.rds")

save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/quality_control/human_cortical_cell_atlas"
seu_type = "cortical"
seu_orig[[]]$tissue_cell = paste0(seu_orig[[]]$Region, ":", seu_orig[[]]$cell_type)
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


### Retain only top 15000 variable genes
obj = brain
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 15000)
saveRDS(obj, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/human_cortical_cell_atlas/MTG_DFC_AnG_orig_comb_hvg.rds")


### Downsample number of cells for each cell type (regardless of tissue)
obj = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/human_cortical_cell_atlas/MTG_DFC_AnG_orig_comb_counts.rds")
N = 20000
Idents(obj) = "cell_type"
obj_sub = subset(x = obj, downsample = N)
saveRDS(obj_sub, paste0("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/human_cortical_cell_atlas/MTG_DFC_AnG_orig_comb_counts_downsamp_cell-type-", N, ".rds"))