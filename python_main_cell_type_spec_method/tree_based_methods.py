import pandas as pd
import numpy as np
import os, argparse, kneed, pickle, logging

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
from sklearn.model_selection import KFold
from sklearn.inspection import permutation_importance
from utils import load_prot_data

warnings.simplefilter("ignore", RuntimeWarning)

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"


def permute_importance(args, prot_spec_final: pd.DataFrame, atlas_smal: pd.DataFrame):
    """
    Calculate permutation importance of features as a negative control
    """
    # Transform the data
    X_df = atlas_smal.copy()
    na_genes = X_df[X_df.isna().any(axis=1)].index.tolist()
    if len(na_genes) > 0:
        logging.error(f"[WARNING] These genes have NAs in gene expression data!! {na_genes}")
        logging.error(f"[WARNING] These genes will be removed in downstream analyses. "
                    "To fix this, please adjust the gene expression data")

        # Report significant proteomic genes among NA genes
        na_sig_genes = prot_spec_final[
            (prot_spec_final["P_value"] < 5e-7) &
            (prot_spec_final["gene"].isin(na_genes))
        ]["gene"].tolist()
        logging.error(f"[WARNING] Among NA genes, these are those with proteomic pval < 5e-7: {na_sig_genes}")
        X_df = X_df[~X_df.index.isin(na_genes)].copy()
        prot_spec_final = prot_spec_final[~prot_spec_final["gene"].isin(na_genes)].copy()

    # Preprocess some data
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[[col, f"log{col}", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    tmp["z_score"] = (2*(tmp[col] > 1) - 1) * tmp["P_value"].apply(lambda x: stats.norm.isf(x / 2))

    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    hr = tmp[args.output_label]

    if args.abs_hr: hr = np.abs(hr)
    sub_atl = sub_atl.loc[tmp.index, :]
    feature_names = sub_atl.columns.tolist()

    # Split the data k-fold
    kf = KFold(n_splits=args.kfold_n, shuffle=True, random_state=42)
    df_l = []
    for i, (train_index, test_index) in enumerate(kf.split(sub_atl)):

        X_train, y_train = sub_atl.iloc[train_index], hr.iloc[train_index]
        X_test, y_test = sub_atl.iloc[test_index], hr.iloc[test_index]

        if args.ztransform_type == 1:
            obj = StandardScaler()
            obj.fit(X_train)
            X_train, X_test = obj.transform(X_train), obj.transform(X_test) 
        elif args.ztransform_type == 2:
            X_train = StandardScaler().fit_transform(X_train.T).T
            X_test = StandardScaler().fit_transform(X_test.T).T
        else:
            X_train = X_train.to_numpy()
            X_test = X_test.to_numpy()

        # Train the model on the train folds
        model = RandomForestRegressor(
            n_estimators=args.num_trees, min_samples_split=args.min_samples_split, min_samples_leaf=args.min_samples_leaf,
            max_samples=args.max_samples, n_jobs=-1
        )
        model.fit(sub_atl, hr)
            
        train_score = model.score(X_train, y_train)
        test_score = model.score(X_test, y_test)
        logging.error(f"Train score without sample weights: {train_score:.3f}, test score without sample weights: {test_score:.3f}")

        # Validate the model
        r = permutation_importance(model, X_test, y_test, n_repeats=args.n_permute_repeat, random_state=0)
        permute_scores = r.importances   # shape = (n_features, n_repeats)

        # Store results
        df = pd.DataFrame({
            "cell_tissue": feature_names*args.n_permute_repeat, 
            "score": permute_scores.flatten(order="F"),
            "fold": [i]*permute_scores.size
        })
        df_l.append(df)

    # Save results
    final_df = pd.concat(df_l)
    final_df.to_csv(f"{args.save_path}/permute_importance_scores.tsv", sep="\t", index=False)

    # Make plots
    plt.figure(figsize=(9, 20))
    sns.boxplot(data=final_df, x="score", y="cell_tissue", hue="fold")
    plt.axvline(0.0, color='black', linestyle='--')
    plt.title(f"Permutscores at kfold={args.kfold_n}, repeats={args.n_permute_repeat}")
    plt.savefig(f"{args.save_path}/permute_importance_scores.png", bbox_inches="tight", dpi=300)


def hyperparam_search(args, prot_spec_final: pd.DataFrame, atlas_smal_merged: pd.DataFrame):
    """
    Hyperparam search
    """
    na_genes = atlas_smal_merged[atlas_smal_merged.isna().any(axis=1)].index.tolist()
    if len(na_genes) > 0:
        logging.error(f"[WARNING] These genes have NAs in gene expression data!! {na_genes}")
        logging.error(f"[WARNING] These genes will be removed in downstream analyses. "
                    "To fix this, please adjust the gene expression data")

        # Report significant proteomic genes among NA genes
        na_sig_genes = prot_spec_final[
            (prot_spec_final["P_value"] < 5e-7) &
            (prot_spec_final["gene"].isin(na_genes))
        ]["gene"].tolist()
        logging.error(f"[WARNING] Among NA genes, these are those with proteomic pval < 5e-7: {na_sig_genes}")
        atlas_smal_merged = atlas_smal_merged[~atlas_smal_merged.index.isin(na_genes)].copy()
        prot_spec_final = prot_spec_final[~prot_spec_final["gene"].isin(na_genes)].copy()

    # Preproces some data
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"
    obj = StandardScaler()

    atlas_smal_subset = atlas_smal_merged
    if (args.ztransform_type == 1): sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
    elif args.ztransform_type == 2: sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy().T).T
    else: sub_atl = atlas_smal_subset.to_numpy()
    sub_atl = pd.DataFrame(sub_atl, columns=atlas_smal_subset.columns, index=atlas_smal_subset.index)

    logging.error(f"2: {sub_atl.shape}")
    logging.error(sub_atl.head())

    tmp = sub_atl.merge(prot_spec_final[[col, f"log{col}", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    tmp["z_score"] = (2*(tmp[col] > 1) - 1) * tmp["P_value"].apply(lambda x: stats.norm.isf(x / 2))
    logging.error(f"3: {tmp.shape}")

    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    
    hr = tmp[args.output_label]
    if args.abs_hr: hr = np.abs(hr)
    sub_atl = sub_atl.loc[tmp.index, :]

    # List of params to search over
    param_grid = {
        "n_estimators": [100, 200, 500],
        "max_features": ["sqrt", "log2", None],
        "max_depth": [None, 10, 20],
        "min_samples_leaf": [1, 4]
    }

    results = []

    for n in param_grid["n_estimators"]:
        for mf in param_grid["max_features"]:
            for md in param_grid["max_depth"]:
                for ms in param_grid["min_samples_leaf"]:
                    
                    model = RandomForestRegressor(
                        n_estimators=n,
                        max_features=mf,
                        max_depth=md,
                        min_samples_leaf=ms,
                        oob_score=True,
                        bootstrap=True,
                        n_jobs=-1,
                        random_state=42
                    )
                    model.fit(sub_atl, hr)
                    results.append(((n, mf, md, ms), model.oob_score_))

    # Best results
    best_params, best_oob = max(results, key=lambda x: x[1])
    args.num_trees = best_params[0]
    args.max_features = best_params[1]
    args.max_depth = best_params[2]
    args.min_samples_leaf = best_params[3]

    # Save the results
    tmp = [[*i[0], i[1]] for i in results]
    df = pd.DataFrame(tmp, columns=["n_estimators", "max_features", "max_depth", "min_sample_leaf", "r2"])
    df.to_csv(f"{args.save_path}/hyperparam_search_results.tsv", sep="\t", index=False)

    return args


def random_forests(args, prot_spec_final: pd.DataFrame, atlas_smal_merged: pd.DataFrame):
    """
    Run random forests
    """
    na_genes = atlas_smal_merged[atlas_smal_merged.isna().any(axis=1)].index.tolist()
    if len(na_genes) > 0:
        logging.error(f"[WARNING] These genes have NAs in gene expression data!! {na_genes}")
        logging.error(f"[WARNING] These genes will be removed in downstream analyses. "
                    "To fix this, please adjust the gene expression data")

        # Report significant proteomic genes among NA genes
        na_sig_genes = prot_spec_final[
            (prot_spec_final["P_value"] < 5e-7) &
            (prot_spec_final["gene"].isin(na_genes))
        ]["gene"].tolist()
        logging.error(f"[WARNING] Among NA genes, these are those with proteomic pval < 5e-7: {na_sig_genes}")
        atlas_smal_merged = atlas_smal_merged[~atlas_smal_merged.index.isin(na_genes)].copy()
        prot_spec_final = prot_spec_final[~prot_spec_final["gene"].isin(na_genes)].copy()

    # Preprocess data
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"
    obj = StandardScaler()

    atlas_smal_subset = atlas_smal_merged
    if (args.ztransform_type == 1): sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
    elif args.ztransform_type == 2: sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy().T).T
    else: sub_atl = atlas_smal_subset.to_numpy()
    sub_atl = pd.DataFrame(sub_atl, columns=atlas_smal_subset.columns, index=atlas_smal_subset.index)

    tmp = sub_atl.merge(prot_spec_final[[col, f"log{col}", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    tmp["z_score"] = (2*(tmp[col] > 1) - 1) * tmp["P_value"].apply(lambda x: stats.norm.isf(x / 2))

    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    
    hr = tmp[args.output_label]
    if args.abs_hr: hr = np.abs(hr)
    sub_atl = sub_atl.loc[tmp.index, :]

    # Random forest
    model = RandomForestRegressor(
        n_estimators=args.num_trees, min_samples_split=args.min_samples_split, min_samples_leaf=args.min_samples_leaf,
        max_samples=args.max_samples, max_features=args.max_features, max_depth=args.max_depth, n_jobs=-1
    )
    model.fit(sub_atl, hr)

    # Save results
    scores = model.feature_importances_  # one-dimensional numpy array
    coef_tree_df = pd.DataFrame({"tree_feature_importance": scores, "cell_tissue": sub_atl.columns}).sort_values(by="tree_feature_importance", ascending=False)
    coef_tree_df.to_csv(f"{args.save_path}/coef_random_forest.tsv", sep="\t", index=False)
    with open(f"{args.save_path}/random_forest_model.pkl", "wb") as f:
        pickle.dump(model, f) 
    

def main(args):

    # Save the arguments
    os.makedirs(args.save_path, exist_ok=True)

    # Load in the atlas data
    full_atlas = pd.read_csv(args.atlas_path, sep="\t")
    if "gene" in full_atlas.columns: full_atlas = full_atlas.set_index("gene")
    atlas_smal = pd.read_csv(args.atlas_smal_path, sep="\t").set_index("gene")

    # Load in prot data
    logging.error("Loading prot data...")
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)
    prot_spec_final.to_csv(f"{args.save_path}/prot_spec_final.tsv", sep="\t", index=False)

    # Convert some argument values to bool
    args.abs_hr = args.abs_hr == 1

    # If param search
    if (args.param_search == 1):
        logging.error("Running hyperparam search...")
        args = hyperparam_search(args, prot_spec_final, atlas_smal)
    
    # Save params
    args_dict = vars(args)
    with open(f"{args.save_path}/cmd_args.json", "w") as json_file:
        json.dump(args_dict, json_file, indent=4)

    # Train
    logging.error("Training the model...")
    random_forests(args, prot_spec_final, atlas_smal)

    # Run permutation importance for random forests (more reliable than impurity-based feature importances)
    logging.error("Permutation importance...")
    permute_importance(args, prot_spec_final, atlas_smal)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas_path", type=str, help="Atlas path")
    parser.add_argument("--atlas_smal_path", type=str, help="Atlas smal path")
    parser.add_argument("--prot_data_path", type=str, help="Prot path")
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--disease", type=str)
    parser.add_argument("--output_label", type=str, default="HR")
    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to set HR to abs value")

    parser.add_argument("--param_search", type=int, default=0, help="Whether to perform param search instead")    
    parser.add_argument("--num_trees", type=int, default=500, help="Number of trees")
    parser.add_argument("--min_samples_split", type=int, default=10, help="Minimum size of a node before getting split")
    parser.add_argument("--min_samples_leaf", type=int, default=2, help="Minimum size of a leaf")
    parser.add_argument("--max_samples", type=float, default=1.0, help="Fraction of dataset for bootstrapping")
    parser.add_argument("--max_depth", type=int, default=None, help="Maximum depth of trees")

    parser.add_argument("--kfold_n", type=int, default=5, help="Number of k for kfolds")
    parser.add_argument("--n_permute_repeat", type=int, default=30, help="Number of n permutations")

    parser.add_argument("--ztransform_type", type=int, default=1, help="Whether to z transform on each cell type (1) or each gene (2)")

    args = parser.parse_args()
    main(args)