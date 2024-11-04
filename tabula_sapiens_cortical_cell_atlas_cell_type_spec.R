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

output_dir = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/results/cell_type_specificity"
atlas = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens_cortical_cells/tabula_sapiens_cortical_cells_atlas_all-genes_downsamp.rds")
atlas = NormalizeData(atlas)

### Cell type specificity by cell
#atlas[[]]$tissue_cell = paste0(atlas[[]]$tissue, ":", atlas[[]]$cell_type)
Idents(atlas) = "cell_type"

num_max_cells = 1000
markers_wilcox_max = FindAllMarkers(
  object = atlas,
  slot = "data",
  test.use = "wilcox",
  assay = "RNA",
  max.cells.per.ident = num_max_cells
)

write.table(
  markers_wilcox_max, 
  file = paste0(output_dir, "/markers_wilcox_cell_max-", num_max_cells, "_downsamp.tsv"), 
  sep = "\t", 
  quote = FALSE
)