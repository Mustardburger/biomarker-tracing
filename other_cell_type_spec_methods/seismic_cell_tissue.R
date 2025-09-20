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
  library("seismicGWAS")
})
options(future.globals.maxSize = 10000 * 1024^2)

# Get the function from seismic, with modifications
get_ct_trait_associations <- function(sscore, magma, magma_gene_col = "GENE",
                                      magma_z_col = "ZSTAT", model = "linear") {
  pvalue <- FDR <- NULL # due to non-standard evaluation notes in R CMD check
  
  # validate method parameter
  if (! model %in% c("linear", "spearman")){
    stop("The model option parameter can be only either 'linear' for linear regression or
         'spearman' for Spearman's rank correlation.")
  }
  
  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- melt(sscore, id.vars = "gene", variable.name = "cell_type", value.name = "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$cell_type <- as.character(sscore$cell_type)
  sscore$specificity <- as.numeric(sscore$specificity)
  
  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)
  
  # get all cell type names
  cts <- unique(sscore$cell_type)
  
  # combine the magma annots with the specificities
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  
  # calculate the association for each cell type
  if (model == "linear"){
    res <- rbindlist(lapply(cts, function(ct) {
      slm <- speedglm::speedlm(dt[cell_type == ct]$zstat ~ dt[cell_type == ct]$specificity) # fast lm
      slm_summ <- summary(slm)$coefficients
      
      # extract the 1-sided hypothesis test p-value
      pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4] / 2 else (1 - slm_summ[2, 4] / 2)
      return(data.table(cell_type = ct, pvalue = pval))
    }))
  }else{
    res <- rbindlist(lapply(cts, function(ct) {
      ct_cortest <- stats::cor.test(dt[cell_type == ct]$specificity, dt[cell_type == ct]$zstat, 
                                    alternative = "greater", method = "spearman", exact = FALSE)
      
      # extract p values and avoid 0
      pval <- if (ct_cortest$p.value > 0) ct_cortest$p.value else 2.2e-16 
      return(data.table(cell_type = ct, pvalue = pval))
    }))
  }
  
  res[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  res <- res[order(pvalue, FDR)]
  
  return(res)
}

# Do some preprocessing
atlas_smal = fread()
prot_spec_final = fread()



# Read in the data
obj = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data/tabula_sapiens/tabula_sapiens_all_cells.rds")
# Idents(obj) = "cell_type"
# print(obj)
# obj = subset(obj, idents = c("unknown"), invert = TRUE)
# #obj = NormalizeData(obj)
# 
# ### seismic cell type specificity score ###
# DefaultAssay(obj) = "RNA"
# 
# obj_diet = DietSeurat(obj, assays = "RNA", layers = "data")
# obj_sce = as.SingleCellExperiment(obj_diet)
# obj_sscore <- calc_specificity(obj_sce, ct_label_col='cell_type', min_avg_exp_ct=0.0)
# 
# file_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/results/cell_type_specificity/tabula_sapien_all_only_seismic.tsv"
# write.table(as.data.frame(obj_sscore), file = file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)