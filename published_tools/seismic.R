#' Scripts to generate seismic cell type specificity matrix

suppressMessages({
  library("plyr")
  library("Seurat")
  library("qs")
  library("dplyr")
  library("tidyr")
  library("data.table")
  # library('seismicGWAS')
  library("stringr")
  library("glue")
})

#' get_ct_trait_associations
#' Rewrite this function a bit to fit with proteomic input
#'
#' @param sscore: Matrix from calc_specificity
#' @param magma: data.table from prot_df
get_ct_trait_associations <- function(sscore, magma, magma_gene_col = "GENE",
                                      magma_z_col = "ZSTAT") {
  pvalue <- FDR <- NULL # due to non-standard evaluation notes in R CMD check

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
  res <- rbindlist(lapply(cts, function(ct) {
    slm <- speedglm::speedlm(dt[cell_type == ct]$zstat ~ dt[cell_type == ct]$specificity) # fast lm
    slm_summ <- summary(slm)$coefficients

    # extract the 1-sided hypothesis test p-value
    pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4] / 2 else (1 - slm_summ[2, 4] / 2)
    return(data.table(cell_type = ct, pvalue = pval))
  }))

  res[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  res <- res[order(pvalue, FDR)]

  return(res)
}


# Paths
seurat_obj_path = ""    # Seurat path
thres = 50              # Minimum number of cells in each cell type
gene_path = ""          # Path to a dataframe with $gene as column, $gene containing all proteomic genes
cell_spec_save_path = ""  # Path to save seismic cell type specificity table
run_name = ""           # Run name

# Read in seurat obj
seu = qread(seurat_obj_path)

# Create cell_tissue column
df = seu[[]]
df <- df %>%
  dplyr::mutate(
    cell_tissue = paste0(
      str_replace_all(tissue, "_", "."), "_", str_replace_all(cell_type, " ", ".")
    )
  )
seu$cell_tissue = df$cell_tissue

# Remove cell_tissue that are smaller than 50 cells
l = table(df$cell_tissue); l = l[l <= thres]; to_exclude = names(l)
seu_smal = subset(seu, subset = cell_tissue %in% to_exclude, invert = TRUE)

# Grab the list of plasma genes
gene_names = fread(gene_path)$gene

# Calculate seismic cell type specificity
seu_smal = NormalizeData(seu_smal)
seu_diet = DietSeurat(seu_smal, features = gene_names)
seu_sce = as.SingleCellExperiment(seu_diet)
seismic_spec_scores <- calc_specificity(
  seu_sce, ct_label_col='cell_tissue', 
  min_avg_exp_ct=0.0, min_cells_gene_exp=0, min_uniq_ct=0,
  min_ct_size = 0
)
seismic_df = as.data.frame(seismic_spec_scores)
seismic_df$gene = rownames(seismic_df)
rownames(seismic_df) = NULL

# Using seismic calc_specificity as is, there are genes that are in the SCE object
# but are all NaN in seismic's specificity score (e.g. ENSG00000243509 for HPA data)
# that is because in the original SCE object the gene has 0 expression across all cells
fwrite(seismic_df, cell_spec_save_path, 
  sep="\t", quote=FALSE, row.names=FALSE
)