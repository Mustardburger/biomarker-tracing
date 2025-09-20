import argparse, logging
import os, json
import numpy as np
import pandas as pd
from scipy.stats import norm, shapiro, t
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from sklearn.preprocessing import StandardScaler
from utils import *

def _gini_coeff(df, g):
    arr = df.loc[g, :].to_numpy()
    if np.all(arr == 0):
        return np.nan

    arr = arr[arr > 0.0]
    
    x = np.sort(arr)
    n = len(x)
    cumx = np.cumsum(x)
    
    gini = (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return gini


def compute_covariates(args, atlas_smal):
    """
    atlas_smal has shape (num genes, num cell tissues)
    with the index being the gene names
    """
    gene_l = None

    # Calculate Gini coefficient
    if args.covar_gini == 1:
        logging.error("Calculating Gini scores:...")
        gene_l = atlas_smal.index.tolist()
        gini_l = [_gini_coeff(atlas_smal, g) for g in gene_l]
    
    if gene_l is not None:
        covar_df = pd.DataFrame({
            "gene": gene_l,
            "gini": gini_l
        }).set_index("gene")
        return covar_df
    else:
        return None


def univariate_testing(args, atlas_smal, prot_spec_final, covar_df=None):
    """
    Make sure atlas_smal, prot_spec_final, and covar_df both have gene as index
    """
    with open(f"{args.save_path}/model_summary.txt", "a") as f:
        f.write("Linear model summary: \n")

    # do some prep
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"
    prot_spec_final, _, _ = prep_data(args, prot_spec_final[[col, f"log{col}", "P_value"]], col)

    # atlas_smal has shape = (num genes, num cell types)
    obj = StandardScaler()
    scaled = obj.fit_transform(atlas_smal)
    result_df = pd.DataFrame(scaled, columns=atlas_smal.columns, index=atlas_smal.index)

    # Because of python string format, replace all . with _
    old_col_names = result_df.columns.tolist()
    new_col_names = [i.replace(".", "_") for i in result_df.columns.tolist()]
    result_df.columns = new_col_names
    cell_tis = new_col_names

    # Add covariates if any
    covar_cols = []
    if covar_df is not None:
        result_df = result_df.merge(covar_df, right_index=True, left_index=True)
        covar_cols = covar_df.columns.tolist()
    prot_df_sub = prot_spec_final[[args.output_label]]
    if args.abs_hr == 1:
        prot_df_sub[args.output_label] = np.abs(prot_df_sub[args.output_label])
    result_df = result_df.merge(prot_df_sub, right_index=True, left_index=True)

    logging.error(covar_cols)
    logging.error(covar_df)

    # Loop over each cell-tissue in the dataset
    final_df = []
    for old_ct, ct in zip(old_col_names, cell_tis):

        indep_var = ct
        all_cols = covar_cols + [indep_var, args.output_label]
        logging.error(f"Processing {ct}")
        result_df_ct = result_df[all_cols].dropna()

        covariates = covar_cols
        formula = args.output_label + " ~ " + " + ".join([indep_var] + covariates)

        # Linear regression
        model = smf.ols(formula, data=result_df_ct).fit()
        with open(f"{args.save_path}/model_summary.txt", "a") as f:
            f.write(f"Linear model for {old_ct}\n")
            f.write(model.summary().as_text())
            f.write("\n===============================\n")

        # Residual normality test
        shapiro_p = shapiro(model.resid)[1]

        # Extract coefficients
        coeffs = model.params
        tvals = model.tvalues
        pvals = model.pvalues

        beta_indep = coeffs[indep_var]
        tval_indep = tvals[indep_var]
        pval_indep = pvals[indep_var]

        intc = coeffs["Intercept"]
        tval_intc = tvals["Intercept"]
        pval_intc = pvals["Intercept"]

        r2 = model.rsquared

        # One-sided test
        df_resid = model.df_resid
        pval_one_sided = 1 - t.cdf(tval_indep, df=df_resid)

        # Append
        final_df.append([old_ct, beta_indep, pval_indep, tval_indep, intc, pval_intc, tval_intc, r2, shapiro_p, pval_one_sided])

    ### Save results
    final_df = pd.DataFrame(final_df, columns=[
        "cell_tissue", "beta_predictor", "pval_predictor", "tval_predictor",
        "intercept", "pval_intercept", "tval_intercept", "r2_score",
        "residual_pval", "pval_predictor_one_side"]
    )

    # FDR correction
    final_df["fdr_predictor"] = multipletests(final_df["pval_predictor"], method="fdr_bh")[1]

    # FDR correction one-side, removing negative coeff samples first
    pos_coeff = final_df[final_df["beta_predictor"] > 0.0].copy()[["cell_tissue", "pval_predictor_one_side"]].copy()
    pos_coeff["fdr_one_side_predictor"] = multipletests(pos_coeff["pval_predictor_one_side"], method="fdr_bh")[1]
    final_df = final_df.merge(pos_coeff.drop(columns="pval_predictor_one_side"), on="cell_tissue", how="left").fillna(1.0).sort_values("fdr_one_side_predictor")

    # Write output
    final_df.to_csv(
        os.path.join(args.save_path, f"univar_regression_results.tsv"),
        sep="\t", index=False
    )
    prot_df_sub.to_csv(
        os.path.join(args.save_path, "prot_spec_final_sub.tsv"),
        sep="\t"
    )


def main(args):

    if args.covar_df == "None": args.covar_df = None

    # Save the arguments
    os.makedirs(args.save_path, exist_ok=True)
    args_dict = vars(args)
    with open(f"{args.save_path}/cmd_args.json", "w") as json_file:
        json.dump(args_dict, json_file, indent=4)

    # Load in the atlas data
    full_atlas = pd.read_csv(args.atlas_path, sep="\t")
    if "gene" in full_atlas.columns: full_atlas = full_atlas.set_index("gene")
    atlas_smal = pd.read_csv(args.atlas_smal_path, sep="\t").set_index("gene")

    # Load in prot data
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)
    if "gene" in prot_spec_final.columns: prot_spec_final = prot_spec_final.set_index("gene")
    prot_spec_final.to_csv(f"{args.save_path}/prot_spec_final.tsv", sep="\t", index=False)

    # Load in covariate data
    covar_df = None
    if args.covar_df is not None: 
        covar_df = pd.read_csv(args.covar_df, sep="\t")
        if "gene" in covar_df.columns: covar_df = covar_df.set_index("gene")
    covar_df_more = compute_covariates(args, atlas_smal)
    if covar_df is None and covar_df_more is None: pass
    else: covar_df = pd.concat([covar_df, covar_df_more], axis=1).fillna(0.0)

    # Convert some argument values to bool
    args.abs_hr = args.abs_hr == 1

    # Run univariate testing
    univariate_testing(args, atlas_smal, prot_spec_final, covar_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atlas_path", type=str, help="Atlas path")
    parser.add_argument("--atlas_smal_path", type=str, help="Atlas smal path")
    parser.add_argument("--prot_data_path", type=str, help="Prot path")
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--disease", type=str)
    parser.add_argument("--output_label", type=str, default="HR")
    parser.add_argument("--abs_hr", type=int, default=0,
                        help="Whether to take the absolute value of the output")
    parser.add_argument("--gene_weight_minmax", type=int, default=0, help="Whether to minmax gene weight")

    parser.add_argument("--covar_df", type=str, default=None)
    parser.add_argument("--covar_gini", type=int, default=0)

    args = parser.parse_args()
    main(args)
