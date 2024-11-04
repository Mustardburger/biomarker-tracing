suppressMessages({
  #library("scran")
  library("plyr")
  library("ggrepel")
  library("ggthemes")
  library("grid")
  library("Seurat")
  library("SeuratDisk")
  library("SeuratData")
  library("scales")
  library("qs")
  library("scRNAseq")
  library("dplyr")
  library("tidyr")
  library("data.table")
  library("ggplot2")
  #library("scCustomize")
  library("zellkonverter")
  library('seismicGWAS')
})

source("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/helper_scripts/gazestani_2023_massive.R")
source("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/helper_scripts/general.R")


brain_piehl = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/piehl_nph_rerun_finer/sct_merged_piehl_nph_brain_subtypes.rds")
brain_piehl = brain_piehl[, brain_piehl@meta.data$cell != "Undetermined"]
feature_names <- rownames(brain_piehl@assays$SCT@counts)


### seismic cell type specificity score ###
# Using common genes
common_genes_txt = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seurat_funcs/common_genes_nph_finer_piehl.txt"
common_genes = readLines(common_genes_txt)
DefaultAssay(brain_piehl) = "SCT"

brain_piehl_sct_diet = DietSeurat(
  brain_piehl, assays = "SCT", features = common_genes
)

brain_piehl_sce = as.SingleCellExperiment(brain_piehl_sct_diet)
# head(colData(brain_piehl_sce))
brain_piehl_sscore <- calc_specificity(brain_piehl_sce, ct_label_col='cell', min_avg_exp_ct=0.0)

file_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic/nph_piehl_brain_subtypes_seismic_common_genes.tsv"
write.table(as.data.frame(brain_piehl_sscore), file = file_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
brain_piehl_sscore = read.table(file_path, sep = "\t", header = TRUE)

### seismic influential gene analysis ###
# Rewrite the function find_inf_genes from seismicGWAS
find_inf_genes_for_proteomics <- function(ct, sscore, magma, metric_col = "AUC.ROC") {

  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- sscore[, c("gene", ct), with = F]
  names(sscore) <- c("gene", "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$specificity <- as.numeric(sscore$specificity)

  magma <- magma[, c("gene", metric_col)]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)

  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  lm_out <- stats::lm(zstat ~ specificity, data = dt)
  lm_summ <- summary(lm_out)$coefficients
  lm_pval <- if (lm_summ[2, 1] > 0) lm_summ[2, 4] / 2 else (1 - lm_summ[2, 4] / 2)
  if (lm_pval < 0.05) {
    cat("Cell type ", ct, " does seem to have a significant association with trait
         (estimated p-value ", round(lm_pval, 4), " based on a 1-sided hypothesis test,
         prior to multiple hypothesis test correction, so influential gene analysis is relevant")
  }

  dt[, dfbetas := stats::dfbetas(lm_out)[, 2]]
  thresh <- 2 / sqrt(nrow(dt))
  # only considering positive dfbetas values (since our analyses is based off of
  # the 1-sided test of positive relationships between specificity + risk)
  dt[, is_influential := ifelse(dfbetas > thresh, T, F)]
  dt <- dt[order(-dfbetas)]

  return(list(df = dt, pval = lm_pval))
}

# Load in the proteomics dataset
rep1_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/avani_work/single_biomarker_performance/data_output/01_single_biomarker/rep1_single_biomarker_auc_phuc_rerun.csv"
emory_MS_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/avani_work/single_biomarker_performance/data_output/02_single_biomarker/emory300MS_phuc_rerun/emory_MS_performance_phuc_rerun.csv"
emory_SOMA_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/avani_work/single_biomarker_performance/data_output/02_single_biomarker/emory300SOMA_phuc_rerun/emory_SOMA_performance_phuc_rerun.csv"

# Some preprocessing
prot_df_path = emory_MS_path
prot_df = read.table(prot_df_path, sep = ",", header = TRUE)
prot_df = subset(prot_df, select = -X)
prot_df <- prot_df %>%
  separate(Biomarker, into = c("gene", "protID"), sep = "\\|", remove = FALSE)

# Handling dup values
gene_counts = table(prot_df$gene)
dup_genes = names(gene_counts[gene_counts > 1])
mean_val_dup <- aggregate(cbind(AUC.ROC, Accuracy) ~ gene, data = prot_df[prot_df$gene %in% dup_genes, ], FUN = mean)
prot_df_fin <- rbind(
    subset(prot_df[!prot_df$gene %in% dup_genes, ], select = -c(Biomarker, protID)), 
    mean_val_dup
)

# Filter out common genes
common_genes = rownames(brain_piehl_sscore)
sc_prot_genes = intersect(common_genes, prot_df_fin$gene)
prot_df_fin = prot_df_fin[prot_df_fin$gene %in% sc_prot_genes, ]
brain_piehl_sscore = brain_piehl_sscore[rownames(brain_piehl_sscore) %in% sc_prot_genes, ]

# Try running influential genes
inf_genes = list()
pval_l = list()
for (ct in colnames(brain_piehl_sscore)) {
  val = find_inf_genes_for_proteomics(ct, brain_piehl_sscore, prot_df)
  inf_genes[[ct]] <- val$df
  pval_l[[ct]] = val$pval
}

# Save results
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic"
save_prefix = "emory_MS"
inf_genes = readRDS(paste0(save_path, "/", save_prefix, "_influential_genes_piehl_nph_brain_subtypes.rds"))
pval_l = readRDS(paste0(save_path, "/", save_prefix, "_seismic_corr_pval_piehl_nph_brain_subtypes.rds"))

# Sort out cell types with pval smaller than 0.05
sig_cells <- names(pval_l)[sapply(pval_l, function(x) x < 0.05)]
cell_names = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/piehl_nph_rerun_finer/cell_names_piehl_nph.tsv", sep = "\t", header = TRUE)
cell_names = cell_names[cell_names$clusters %in% sig_cells, ]
full_names = paste0(cell_names$cell, "_", cell_names$clusters)
piehl_cells = setdiff(sig_cells, cell_names$clusters)
sig_cells = c(piehl_cells, full_names)

cell = "PVALB_HOMER3"
a = inf_genes[[cell]][inf_genes[[cell]]$is_influential == TRUE, ]
paste(sort(a$gene), collapse = ",")

inf_genes[["PVALB_NOG"]][inf_genes[["PVALB_NOG"]]$is_influential == TRUE, ] # ToppGene results interesting
inf_genes[["NDNF_TGFBR2"]][inf_genes[["NDNF_TGFBR2"]]$is_influential == TRUE, ] # ToppGene results interesting
inf_genes[["NDNF_ARHGAP36"]][inf_genes[["NDNF_ARHGAP36"]]$is_influential == TRUE, ]   # ToppGene results interesting

# Save results
saveRDS(inf_genes, paste0(save_path, "/", save_prefix, "_influential_genes_piehl_nph_brain_subtypes.rds"))
saveRDS(pval_l, paste0(save_path, "/", save_prefix, "_seismic_corr_pval_piehl_nph_brain_subtypes.rds"))



### Compare gene markers produced by seismic and wilcoxon
seismic_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic/nph_piehl_brain_subtypes_seismic_common_genes.tsv"
wilcox_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seurat_funcs/markers_wilcox_nomax_brain_subtypes.tsv"
seismic = read.table(file = seismic_path, sep = "\t")
wilcox = read.table(file = wilcox_path, sep = "\t", header = TRUE)

# Rename a bunch of faulty cell types
names(seismic)[names(seismic) == "B.Cells"] <- "B Cells"
names(seismic)[names(seismic) == "CD14..CD16..CD68hi.Monoctyes"] <- "CD14+/CD16+/CD68hi Monoctyes"
names(seismic)[names(seismic) == "CD14..CD16..CD68mid.Monocytes"] <- "CD14+/CD16+/CD68mid Monocytes"
names(seismic)[names(seismic) == "CD14..CD16..CD68lo.Monocytes"] <- "CD14+/CD16-/CD68lo Monocytes"
names(seismic)[names(seismic) == "CD4..T.Cells"] <- "CD4+ T Cells"
names(seismic)[names(seismic) == "CD4..CD8..T.Cells"] <- "CD4+/CD8+ T Cells"
names(seismic)[names(seismic) == "CD8..T.Cells"] <- "CD8+ T Cells"
names(seismic)[names(seismic) == "T.Regulatory.Cells"] <- "T Regulatory Cells"
names(seismic)[names(seismic) == "NK.Cells"] <- "NK Cells"
names(seismic)[names(seismic) == "FEZF2_RP11.360K13.1"] <- "FEZF2_RP11-360K13.1"

cell_types = unique(wilcox$cluster)
seismic_elimi = 1e-7
wilcox_elimi = 0.0
wilcox_pval_thres = 0.01
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seurat_wilcox_seismic"

for (cell_ex in cell_types) {

  if (!(cell_ex %in% colnames(seismic))) {
    print(paste0(cell_ex, " not in seismic!!"))
    next
  }

  wilcox_ex = wilcox[wilcox$cluster == cell_ex, c("p_val", "avg_log2FC", "gene")]
  seismic_ex = seismic[, c(cell_ex), drop = FALSE]
  names(seismic_ex)[names(seismic_ex) == cell_ex] <- "seismic"
  seismic_ex$gene = rownames(seismic_ex)

  wil_sei = merge(wilcox_ex, seismic_ex, by = "gene", all.x = TRUE)
  wil_sei <- na.omit(wil_sei)
  wil_sei = subset(wil_sei, seismic > seismic_elimi & avg_log2FC > wilcox_elimi)
  wil_sei$wilcox_sig = wil_sei$p_val < wilcox_pval_thres

  model <- lm(seismic ~ avg_log2FC, data = wil_sei)
  summary_model <- summary(model)
  r_squared <- summary_model$r.squared

  g = ggplot(wil_sei, aes(x = avg_log2FC, y = seismic, color = wilcox_sig)) +
    geom_point(color = "blue") +
    labs(title = paste0("wilcox avg log2FC vs seismic score in ", cell_ex, "\n r2 = ", round(r_squared, 3)), x = "avg_log2FC", y = "seismic") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5), # Adjust title size and align to the middle
      axis.title = element_text(size = 15),              # Adjust axis titles size
      axis.text = element_text(size = 12)                # Adjust axis text size
    )
  print(g)
  ggplot2::ggsave(
    filename = paste0("seismic_wilcox_", cell_ex, ".png"), 
    g, 
    path = save_path, bg = "white"
  )

}


### Compare influential genes found by seissmic and DEGs in single cell objects
# loads for about 10min
nph_deg_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/nph_degs/author_orig_NPH_cellType_DEres.csv"
nph_deg_df = read.table(nph_deg_path, header = TRUE, sep = ",")

# if the above takes too long, then try this:
nph_deg_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/gazestani_2023_sc_massive/NPH_DE_res/NPH_cellType_DEres.qs"
nph_deg_df = qread(nph_deg_path)

# seismic data
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic"
save_prefix = "emory_MS"
inf_genes = readRDS(paste0(save_path, "/", save_prefix, "_influential_genes_piehl_nph_brain_subtypes.rds"))
pval_l = readRDS(paste0(save_path, "/", save_prefix, "_seismic_corr_pval_piehl_nph_brain_subtypes.rds"))

# Sort out cell types with pval smaller than 0.05
sig_cells <- names(pval_l)[sapply(pval_l, function(x) x < 0.05)]
cell_names = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/piehl_nph_rerun_finer/cell_names_piehl_nph.tsv", sep = "\t", header = TRUE)
cell_names = cell_names[cell_names$clusters %in% sig_cells, ]
full_names = paste0(cell_names$cell, "_", cell_names$clusters)
piehl_cells = setdiff(sig_cells, cell_names$clusters)
sig_cells = c(piehl_cells, full_names)

# In the publication and [Terra notebook](https://app.terra.bio/#workspaces/fbrihuman/sconline_integrative_analysis/analysis/launch/NPH_AD_differential_expression.ipynb), 
# the authors seem to define dysregulated genes as those with the column `adj.P.Val` < 0.1. 
# The relevant function `.integ.DEvis` is defined in `integ_source_code.R` in the author's codebase.
nph_deg_sig = nph_deg_df[nph_deg_df$adj.P.Val < 0.1, ]

# loop through the cell names and extract common genes
common_genes_l = list()
for (cell in full_names) {
  print(cell)
  nph_deg_genes = unique(nph_deg_sig[nph_deg_sig$name == cell, ]$gene)
  ctype = nph_deg_sig[nph_deg_sig$name == cell, ]$ctype[1]
  seismic_influ_genes = unique(inf_genes[[ctype]][inf_genes[[ctype]]$is_influential == TRUE, ]$gene)
  common_genes = intersect(nph_deg_genes, seismic_influ_genes)
  common_genes_l[[cell]] = common_genes
}

save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic"
save_prefix = "emory_MS"
saveRDS(common_genes_l, paste0(save_path, "/", save_prefix, "_seismic_influ_nph_deg_common_genes.rds"))
common_genes_l = readRDS(paste0(save_path, "/", save_prefix, "_seismic_influ_nph_deg_common_genes.rds"))

# Loop through each cell type in common_genes_l and make scatterplot
nph_deg_col = "jk_logFC_median"
categories = c("AbetaTau", "Abeta 1", "Abeta 2+3")

library("cowplot")
library("ggpubr")
for (cell in names(common_genes_l)) {
  if (length(common_genes_l[[cell]]) == 0) {
    next
  }

  cell_n = sub("^[^_]*_", "", cell)
  seismic = inf_genes[[cell_n]]
  seismic = seismic[seismic$gene %in% common_genes_l[[cell]], ][, c("gene", "dfbetas")]

  # Since the NPH DEG df has many columns, let's pick one column that seems most reasonable and move on from there
  plots = list()
  i = 0
  for (category in categories){
    i = i+1
    nph_deg_l = nph_deg_df[(nph_deg_df$gene %in% common_genes_l[[cell]]) & (nph_deg_df$name == cell) & (nph_deg_df$category == category), ][, c("gene", nph_deg_col, "category", "name")]
    merged = merge(seismic, nph_deg_l)

    model <- lm(as.formula(paste0("dfbetas ~ ", nph_deg_col)), data = merged)
    summary_model <- summary(model)
    r_squared <- summary_model$r.squared
    pval = summary_model$coefficients[, 4][[1]]

    p = ggplot(merged, aes_string(x = "dfbetas", y = nph_deg_col)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    #stat_cor(aes(label = ..r.label..), label.x = 3, label.y = 30) +
    ggtitle(paste0(category, " r2=", sprintf("%.3f", r_squared), " pval=", sprintf("%.8f", pval)))

    plots[[i]] = p
  }
  combined_plot <- plot_grid(plotlist = plots, nrow = 1)
  ggplot2::ggsave(
    filename = paste0(cell, "-", nph_deg_col, ".png"),
    combined_plot,
    width = 20, height = 5,
    path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic/seismic_influ_genes_vs_nph_deg", bg = "white"
  )
}


### Compare influential genes found by seismic and other proteomics papers
# Paper: Blood protein assessment of leading incident diseases and mortality in the UK Biobank
# Supplementary Table 13
uk_bio_prot_df_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/nature_proteomics_papers/blood_protein_assessment_of_leading_incident_diseases_and_mortality_in_the_UK_Biobank/supp_table_13.tsv"
uk_bio_prot_df <- read.delim(uk_bio_prot_df_path, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, fill = TRUE)


# seismic data
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic"
save_prefix = "emory_MS"
inf_genes = readRDS(paste0(save_path, "/", save_prefix, "_influential_genes_piehl_nph_brain_subtypes.rds"))
pval_l = readRDS(paste0(save_path, "/", save_prefix, "_seismic_corr_pval_piehl_nph_brain_subtypes.rds"))
sig_cells <- names(pval_l)[sapply(pval_l, function(x) x < 0.05)]

# check similar genes between uk biobank and seismic influential genes
uk_bio_ad = uk_bio_prot_df[uk_bio_prot_df$ProteinScore == "Alzheimer's dementia", ]
uk_bio_ad_genes = strsplit(uk_bio_ad$Proteins, ", ")[[1]]
for (cell in sig_cells) {
  inf_genes_cell = inf_genes[[cell]][inf_genes[[cell]]$is_influential == TRUE, ]$gene
  number_intersect = length(intersect(uk_bio_ad_genes, inf_genes_cell))
  print(paste0("Number of shared genes between seismic's ", cell, " and UK Biobank AD: ", number_intersect))
}

# apparently no genes showed up. because this list of UK Biobank genes are not present in original gene set
all_seismic_genes = inf_genes[["FEZF2_PDLIM1"]]$gene
number_intersect = length(intersect(uk_bio_ad_genes, all_seismic_genes))
print(paste0("Number of shared genes between all seismic genes and UK Biobank AD: ", number_intersect))
paste(sort(intersect(uk_bio_ad_genes, all_seismic_genes)), collapse = ",")

# since there's no common gene to AD, let's look at other diseases
diseases = unique(uk_bio_prot_df$ProteinScore)
for (disease in diseases) {
  uk_bio_dis = uk_bio_prot_df[uk_bio_prot_df$ProteinScore == disease, ]
  uk_bio_dis_genes = strsplit(uk_bio_dis$Proteins, ", ")[[1]]
  number_intersect = length(intersect(uk_bio_dis_genes, all_seismic_genes))
  print(paste0("Number of genes in UK Biobank ", disease, ": ", length(uk_bio_dis_genes)))
  print(paste0("Number of shared genes between all seismic genes and UK Biobank ", disease, ": ", number_intersect))  
}


### Associations of different gene sets via cell type specificity

# since common genes between gazestani-piehl and other gene sets are so few
# rerun seismic on only gazestani dataset
brain_piehl = readRDS("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/piehl_nph_rerun_finer/sct_merged_piehl_nph_brain_subtypes.rds")
brain_piehl = brain_piehl[, brain_piehl@meta.data$dataset == "brain"]
DefaultAssay(brain_piehl) = "SCT"
brain_piehl_sct_diet = DietSeurat(brain_piehl, assays = "SCT")
brain_piehl_sce = as.SingleCellExperiment(brain_piehl_sct_diet)
brain_piehl_sscore <- calc_specificity(brain_piehl_sce, ct_label_col='cell', min_avg_exp_ct=0.0)
write.table(brain_piehl_sscore, "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic/nph_finer_seismic_all_genes.tsv", sep = "\t", quote = FALSE)
brain_piehl_sscore = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic/nph_finer_seismic_all_genes.tsv", sep = "\t", header = TRUE)

# seismic cell type specificity score
seismic_spec_df_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic/nph_finer_piehl_seismic_all_genes.tsv"
seismic_spec_df = read.table(seismic_spec_df_path, sep = "\t", header = TRUE)

# gene list from nature proteomics
uk_bio_prot_df_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/nature_proteomics_papers/blood_protein_assessment_of_leading_incident_diseases_and_mortality_in_the_UK_Biobank/supp_table_13.tsv"
uk_bio_prot_df <- read.delim(uk_bio_prot_df_path, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, fill = TRUE)
uk_bio_ad = uk_bio_prot_df[uk_bio_prot_df$ProteinScore == "Alzheimer's dementia", ]
uk_bio_ad_genes = strsplit(uk_bio_ad$Proteins, ", ")[[1]]

# gene list from seismic analysis
save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/cell_type_specificity/seismic"
save_prefix = "emory_MS"
inf_genes = readRDS(paste0(save_path, "/", save_prefix, "_influential_genes_piehl_nph_brain_subtypes.rds"))
pval_l = readRDS(paste0(save_path, "/", save_prefix, "_seismic_corr_pval_piehl_nph_brain_subtypes.rds"))
sig_cells <- names(pval_l)[sapply(pval_l, function(x) x < 0.05)]

# gene list from csf proteomics paper
csf_prots_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/biomarkers_gene_id.tsv"
csf_prots_df = read.table(csf_prots_path, sep = "\t", header = TRUE)
length(intersect(csf_prots_df$hgnc_symbol, rownames(seismic_spec_df)))

# gene list from seismic paper AD case study
seismic_clin_mg = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/seismic_supplemental_materials/seismic_clin_microglia.tsv", sep = "\t", header = TRUE)
seismic_clin_mg = unique(seismic_clin_mg[seismic_clin_mg$influential == TRUE, ]$gene.symbol)
seismic_tau_ec1 = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/seismic_supplemental_materials/seismic_tau_ECI.tsv", sep = "\t", header = TRUE)
seismic_tau_ec1 = unique(seismic_tau_ec1[seismic_tau_ec1$influential == TRUE, ]$gene.symbol)
seismic_tau_ec2 = read.table("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/seismic_supplemental_materials/seismic_tau_ECII.tsv", sep = "\t", header = TRUE)
seismic_tau_ec2 = unique(seismic_tau_ec2[seismic_tau_ec2$influential == TRUE, ]$gene.symbol)

## possible ideas for analyses

# correlation matrix
all_genes = intersect(union(csf_prots_df$hgnc_symbol, uk_bio_ad_genes), rownames(brain_piehl_sscore))
csf_annot <- data.frame(matrix("CSF", nrow = length(csf_prots_df$hgnc_symbol), ncol = 1))
rownames(csf_annot) = csf_prots_df$hgnc_symbol
colnames(csf_annot) = "Category"

uk_bio_annot <- data.frame(matrix("UK Biobank", nrow = length(uk_bio_ad_genes), ncol = 1))
rownames(uk_bio_annot) = uk_bio_ad_genes
colnames(uk_bio_annot) = "Category"

seis_mg <- data.frame(matrix("Seismic MG", nrow = length(seismic_clin_mg), ncol = 1))
rownames(seis_mg) = seismic_clin_mg
colnames(seis_mg) = "Category"

seis_tau_ec1 <- data.frame(matrix("Seismic Tau EC1", nrow = length(seismic_tau_ec1), ncol = 1))
rownames(seis_tau_ec1) = seismic_tau_ec1
colnames(seis_tau_ec1) = "Category"

seis_tau_ec2 <- data.frame(matrix("Seismic Tau EC2", nrow = length(seismic_tau_ec2), ncol = 1))
rownames(seis_tau_ec2) = seismic_tau_ec2
colnames(seis_tau_ec2) = "Category"

all_cat_annot = rbind(csf_annot, uk_bio_annot, seis_mg, seis_tau_ec1, seis_tau_ec2)
all_genes = rownames(all_cat_annot)

seismic_all = brain_piehl_sscore[rownames(brain_piehl_sscore) %in% all_genes, ]
corr = cor(t(seismic_all))

library(pheatmap)
my_palette <- colorRampPalette(c("blue", "white", "red"))(50)
png("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/gene_set_relationships/csf_uk-bio_seismic-supp_corr_nph_seismic_cell-spec-score.png", width = 1800, height = 1800)
pheatmap(corr, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         fontsize_row = 2, 
         fontsize_col = 2,
         color = my_palette,
         annotation_row = all_cat_annot,
         breaks = seq(-1, 1, length.out = 51))
dev.off()


# clustermap
genes_of_interest = intersect(union(csf_prots_df$hgnc_symbol, uk_bio_ad_genes), rownames(brain_piehl_sscore))
df_of_interest = seismic_all[rownames(seismic_all) %in% genes_of_interest, ]

my_palette <- colorRampPalette(c("white", "red"))(50)
png("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/gene_set_relationships/csf_uk-bio_heatmap_nph_seismic_cell-spec-score.png", width = 1500, height = 1500, res = 250)
pheatmap(df_of_interest, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         fontsize_row = 3, 
         fontsize_col = 3,
         annotation_row = all_cat_annot,
         color = my_palette,
         breaks = seq(0, max(unlist(df_of_interest)), length.out = 51))
dev.off()

# pca
pca_result <- prcomp(seismic_all, scale. = TRUE)
pca_values <- pca_result$x[, 1:2]
pca_df <- as.data.frame(pca_values)
pca_df$RowNames <- rownames(pca_df)
all_cat_annot_sub = all_cat_annot[rownames(all_cat_annot) %in% rownames(pca_df), , drop = FALSE]
pca_df = merge(pca_df, all_cat_annot_sub, by = "row.names", all = TRUE)

# Create a scatter plot of the PCA values with color based on category
png("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/gene_set_relationships/csf_uk-bio_seismic-supp_pca_nph_seismic_cell-spec-score.png", width = 800, height = 800)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Category, label = Row.names)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA Scatter Plot", x = "PC1", y = "PC2") +
  theme_minimal()
dev.off()



### Associations of different gene sets via Gazestani DEGs

# Let's try first between HC and Abeta 1
category = "Abeta 2+3"
nph_deg_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/gazestani_2023_sc_massive/NPH_DE_res/NPH_cellType_DEres.qs"
nph_deg_df = qread(nph_deg_path)
nph_deg_sig = nph_deg_df[(nph_deg_df$adj.P.Val < 0.1) & (nph_deg_df$category == category), ]
genes_of_interest = intersect(union(csf_prots_df$hgnc_symbol, uk_bio_ad_genes), unique(nph_deg_sig$gene))

nph_deg_sig_interest = nph_deg_sig[nph_deg_sig$gene %in% genes_of_interest, c("name", "gene", "jk_logFC_median")]
nph_deg_sig_w <- as.data.frame(pivot_wider(nph_deg_sig_interest, names_from = name, values_from = jk_logFC_median))
nph_deg_sig_w[is.na(nph_deg_sig_w)] <- 0
rownames(nph_deg_sig_w) = nph_deg_sig_w$gene
nph_deg_sig_w = nph_deg_sig_w[, !names(nph_deg_sig_w) %in% c("gene")]

my_palette <- colorRampPalette(c("blue", "white", "red"))(50)
png(
  paste0("/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/gene_set_relationships/csf_uk-bio_heatmap_nph_deg_logfc_", gsub(" ", "-", category), ".png"), 
  width = 1500, height = 800, res = 250
)
pheatmap(nph_deg_sig_w, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         fontsize_row = 3, 
         fontsize_col = 3,
         annotation_row = all_cat_annot,
         color = my_palette,
         breaks = seq(-2, 2, length.out = 51),
         main = category
         )
dev.off()



