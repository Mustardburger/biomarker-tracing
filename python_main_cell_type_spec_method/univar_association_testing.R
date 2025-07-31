library("dplyr")
library("readr")
library("glue")
library("optparse")
library("stats")

### Input arguments

# Define options
option_list <- list(
  make_option(c("--exp_df"), type="character", default=NULL,
              help="Path to gene expression file"),
  make_option(c("--pc_df"), type="character", default=NULL,
              help="Path to principal component file"),
  make_option(c("--prot_df"), type="character", default=NULL,
              help="Path to proteomic result"),
  make_option(c("--num_pc"), type="integer", default=0,
              help="Number of principal components used in model"),
  make_option(c("--mean_covar"), type="character", default="FALSE",
              help="Whether to use mean exp as a covariance"),
  make_option(c("--save_path"), type="character", default=NULL,
              help="Save path"),

  make_option(c("--dep_var"), type="character", default="z",
              help="Dependent variable"),
  make_option(c("--z_transform"), type="integer", default=0,
              help="Whether to transform the data by z-scoring first")
)

# Parse options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

### Some input variables
n_PCs = opt$num_pc  # The number of principal components included
mean_exp_covar = opt$mean_covar  # Whether to use mean expression as a covariate
dep_var = opt$dep_var; z_transform = opt$z_transform
print(glue("Number of PCs used: {n_PCs}, whether to include mean expression term: {mean_exp_covar}"))
dir.create(opt$save_path, showWarnings = FALSE)

### Read in preprocessed tsv file of gene expression
exp_df = read_tsv(opt$exp_df)

# If z transform
if (z_transform == 1) {
    cell_names = colnames(exp_df)[colnames(exp_df) != "gene"]
    z_score_df = exp_df %>% select(-gene) %>% scale() %>% as.data.frame()
    colnames(z_score_df) = cell_names
    z_score_df$gene = exp_df$gene
    exp_df = z_score_df
}

# Calculate mean gene expression for each gene
mean_exp_df = data.frame(gene=exp_df %>% select(gene), mean_exp=rowMeans(exp_df %>% select(-gene)))

### Read in principal components
pc_df = read_tsv(opt$pc_df)
pc_df = pc_df[, c(1:(n_PCs+1))]

### Read in the proteome results
prot_df = read_tsv(opt$prot_df)

prot_df$beta = log10(prot_df$HR)
prot_df$neg_log10_pval = -log10(prot_df$P_value)
log_half_p = -(prot_df$neg_log10_pval) - log(2)
prot_df$z = sign(prot_df$beta) * qnorm(log_half_p, lower.tail = FALSE, log.p = TRUE)
prot_df = prot_df[, c(dep_var, "gene")]

### Setup for linear regression
result_df = merge(merge(exp_df, pc_df, by="gene"), prot_df, by="gene")
if (mean_exp_covar == "TRUE") {
    result_df = merge(result_df, mean_exp_df, by="gene")
}
cell_tis = colnames(exp_df)

# Placeholder variables
beta_indep_l = c()
pval_indep_l = c()
tval_indep_l = c()
beta_intc_l = c()
pval_intc_l = c()
tval_intc_l = c()
r2_l = c()
new_cell_tis = c()
res_pval_l = c()

for (ct in cell_tis) {
    if (ct == "gene") {next}

    new_cell_tis = c(new_cell_tis, ct)
    # Extract only the relevant columns (gene exp of ct, PCs, HRlog10 
    if (n_PCs > 0) {
        pc_cols = paste0("PC", seq(1, n_PCs))
        all_cols = c(ct, pc_cols, dep_var)
    }
    else {
        all_cols = c(ct, dep_var)
    }

    # If include mean expression in model
    if (mean_exp_covar == "TRUE") {
        all_cols = c(all_cols, "mean_exp")
    }
    result_df_ct = result_df[, all_cols]

    # Set up linear reg formula
    indep_var = ct
    covariates <- setdiff(names(result_df_ct), c(dep_var, indep_var))
    formula <- as.formula(paste(dep_var, "~", paste(c(indep_var, covariates), collapse = " + ")))

    ### Linear regression
    model <- lm(formula, data = result_df_ct)
    summary_model <- summary(model)

    # Coefficients table
    coeffs <- summary_model$coefficients

    # Test of normality of residuals
    shapiro_test <- shapiro.test(residuals(model))
    res_pval = shapiro_test$p.value
    res_pval_l = c(res_pval_l, res_pval)
    
    # Beta coefficient and p-value of independent variable
    beta_indep <- coeffs[ct, "Estimate"]
    pval_indep <- coeffs[ct, "Pr(>|t|)"]
    tval_indep <- coeffs[ct, "t value"]
    
    # Intercept stats
    intc <- coeffs["(Intercept)", "Estimate"]
    pval_intc <- coeffs["(Intercept)", "Pr(>|t|)"]
    tval_intc <- coeffs["(Intercept)", "t value"]
    
    # R-squared
    r2 <- summary_model$r.squared

    # Print results
    res = paste0(ct, " has tval: ", formatC(tval_indep, digits=3, format="f"), " and pval: ", formatC(pval_indep, digits=8, format="f"))
    print(res)

    beta_indep_l = c(beta_indep_l, beta_indep)
    pval_indep_l = c(pval_indep_l, pval_indep)
    tval_indep_l = c(tval_indep_l, tval_indep)
    beta_intc_l = c(beta_intc_l, intc)
    pval_intc_l = c(pval_intc_l, pval_intc)
    tval_intc_l = c(tval_intc_l, tval_intc)
    r2_l = c(r2_l, r2)
}

### Save result
final_df = data.frame(
    cell_tissue = as.character(new_cell_tis),
    beta_predictor = beta_indep_l,
    pval_predictor = pval_indep_l,
    tval_predictor = tval_indep_l,
    intercept = beta_intc_l,
    pval_intercept = pval_intc_l,
    tval_intercept = tval_intc_l,
    r2_score = r2_l,
    residual_pval = res_pval_l
)
final_df$fdr_predictor = p.adjust(final_df$pval_predictor)
final_df = final_df %>%
    dplyr::arrange(fdr_predictor)

head(final_df)

write.table(
    final_df, 
    paste0(opt$save_path, "/", "univar_regression_nPC", n_PCs, "_mean-covar-", mean_exp_covar, ".tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
    prot_df, 
    paste0(opt$save_path, "/", "prot_df.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

