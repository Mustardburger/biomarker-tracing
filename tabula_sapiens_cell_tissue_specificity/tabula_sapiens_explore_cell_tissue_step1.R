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


### Step 1: Trimming certain cell-tissue pairs having count below a threshold


# Read in the data
obj = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells.rds")
Idents(obj) = "cell_type"
print(obj)
obj = subset(obj, idents = c("unknown"), invert = TRUE)
#obj = NormalizeData(obj)


### Check cell_tissue relationship - An example
# Extract counts about cell-tissue pairs
meta = obj[[]]
cell_type_count_by_tissue <- meta %>%
  group_by(cell_type, tissue) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(cell_type, desc(count)) %>%
  mutate(cell_type_tissue = paste0(cell_type, "_", tissue))

### Trimming the original Seurat object
# Calculating statistics with a certain threshold
count_thres = 50
cel_tis_below_thres = cell_type_count_by_tissue %>% filter(count < count_thres)
cel_tis_pair = dim(cell_type_count_by_tissue)[1]
cel_tis_pair_tr = dim(cel_tis_below_thres)[1]
tot_cell_counts_remain = cell_type_count_by_tissue %>% filter(count >= count_thres) %>% summarise(count=sum(count)) %>% .[["count"]]
print(paste0("#unique cell-tissue pairs initially: ", cel_tis_pair))
print(paste0("#removed cell-tissue pairs after trimming: ", cel_tis_pair_tr))
print(paste0("#single cells when using thres = ", count_thres, ": ", tot_cell_counts_remain))

# See what cell types or tissues that get removed
cel_tis_remv = unique(cel_tis_below_thres$cell_type_tissue)
cel_tis_remain = setdiff(unique(cell_type_count_by_tissue$cell_type_tissue), cel_tis_remv)
cell_remain = sapply(cel_tis_remain, function(s) strsplit(s, "_")[[1]][1])
tiss_remain = sapply(cel_tis_remain, function(s) strsplit(s, "_")[[1]][2])
cell_removed = setdiff(unique(cell_type_count_by_tissue$cell_type), cell_remain)
tiss_removed = setdiff(unique(cell_type_count_by_tissue$tissue), tiss_remain)
cell_removed
tiss_removed

# Manual check
cell_type_count_by_tissue %>% filter(cell_type == "ciliated epithelial cell")

# Trim the Seurat object
obj[[]] = obj[[]] %>% mutate(cell_type_tissue = paste0(cell_type, "_", tissue))
Idents(obj) = "cell_type_tissue"
obj_tr = subset(obj, idents = cel_tis_remain)

path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens"
saveRDS(obj_tr, paste0(path, "/tabula_sapiens_", count_thres, "_above_cells_per_cell_tissue.rds"))