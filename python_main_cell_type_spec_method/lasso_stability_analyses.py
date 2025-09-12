import pandas as pd
import numpy as np
import os, argparse, kneed, pickle, logging

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
from sklearn.model_selection import KFold
from utils import *

warnings.simplefilter("ignore", RuntimeWarning)

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"

# A function to clean some code
def get_split_data(df, inds, output_label, abs_hr=False):
    X_train = df.iloc[inds, :]
    hr = X_train[output_label]
    if abs_hr: hr = np.abs(hr)
    return X_train, hr


# Using kneedle to find optimal alpha
def optim_alpha_kneedle(args, perf_df: pd.DataFrame):

    Ss = range(20)
    knee_alphas = []
    for S in Ss:
        kneedle = kneed.KneeLocator(perf_df["alpha"], perf_df["zero_coef_perc"], S=S, curve="concave", direction="increasing")
        knee_alpha = kneedle.knee
        knee_alphas.append(knee_alpha)

    coef_perc_at_alpha = perf_df[perf_df["alpha"].isin(knee_alphas)].copy()
    final_kneedle = kneed.KneeLocator(coef_perc_at_alpha["alpha"], coef_perc_at_alpha["zero_coef_perc"], S=1.0, curve="concave", direction="increasing")
    final_alpha = final_kneedle.knee
    if final_alpha is None:
        print("Find alpha manually")
        final_alpha = coef_perc_at_alpha.drop_duplicates(subset="alpha")["alpha"].median()

    perf_df_tmp = perf_df.copy()
    perf_df_tmp["alpha_dist"] = np.abs(perf_df_tmp["alpha"] - final_alpha)
    optim_l1_ratio = perf_df_tmp.sort_values(by="alpha_dist").head(1)["l1_ratio"].item()

    plt.figure()
    plt.plot(perf_df["alpha"], perf_df["zero_coef_perc"], "bo-")
    plt.plot(coef_perc_at_alpha["alpha"], coef_perc_at_alpha["zero_coef_perc"], "ro-")
    plt.axvline(x=final_alpha, color='g', linestyle='--', label=f'Knee point at alpha={final_alpha:.3f}')
    plt.xlabel("alpha")
    plt.ylabel("Percentage of zero coefficients")
    plt.title("Regularization path")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{args.save_path}/kneedle_alpha.png", dpi=300, bbox_inches="tight")

    return final_alpha, optim_l1_ratio


def hyperparam_search(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):

    # Some evidence that this might work. Let's build this pipeline for all diseases
    # diseases = prot_df["Outcome"].unique().tolist()

    alphas, l1_ratios, scores, coeffs, models, conds, pearson_rs = [], [], [], [], [], [], []
    alphas_l = np.logspace(-1, 0, args.num_alphas)
    l1_ratios_l = [1.0]

    obj = StandardScaler()
    X = obj.fit_transform(atlas_smal_merged.to_numpy())
    X_df = pd.DataFrame(X, columns=atlas_smal_merged.columns, index=atlas_smal_merged.index)

    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[[col, "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    tmp["z_score"] = (2*(tmp[col] > 1) - 1) * tmp["P_value"].apply(lambda x: stats.norm.isf(x / 2))
    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    hr = tmp[args.output_label]
    if args.abs_hr: hr = np.abs(hr)
    sub_atl = sub_atl.loc[tmp.index, :]

    if args.gene_weight_minmax: weight_col = "-log10(pval)_minmax"
    else: weight_col = "-log10(pval)"

    logging.error("Starting hyperparam search...")
    for alpha in alphas_l:
        for l1_ratio in l1_ratios_l:

            # Adjust the alphas carefully, because with full cell-tissue dataset, some alphas do not reach convergence
            model = ElasticNet(l1_ratio=l1_ratio, alpha=alpha, positive=args.positive, fit_intercept=args.intercept, max_iter=5000)

            # Catch whether model converges
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", ConvergenceWarning)
                if args.gene_weight: model.fit(sub_atl, hr, sample_weight=tmp[weight_col].tolist())
                else: model.fit(sub_atl, hr)
                if any(issubclass(warning.category, ConvergenceWarning) for warning in w):
                    print(f">>> {args.disease} with alpha={alpha:.3f} and l1_ratio={l1_ratio:.3f} did not converge")
                    continue
    
            # Score the model
            #pred = model.predict(sub_atl)
            #score = r2_score(hr, pred)
            if args.gene_weight: score = model.score(sub_atl, hr, sample_weight=tmp[weight_col].tolist())
            else: score = model.score(sub_atl, hr)
            pred = model.predict(sub_atl)
            r, _ = pearsonr(hr, pred)
            r = np.nan_to_num(r)
        
            alphas.append(alpha)
            scores.append(score)
            pearson_rs.append(r)
            coeffs.append(model.coef_)
            l1_ratios.append(l1_ratio)
            conds.append(args.disease)
            models.append(model)

    perf_df = pd.DataFrame({"disease": conds, "alpha": alphas, "l1_ratio": l1_ratios, "r2": scores, "pearson_r": pearson_rs})
    coef_np = np.array(coeffs)
    coef_df = pd.DataFrame(coef_np, columns=sub_atl.columns)
    coef_df["disease"] = conds
    perf_df["zero_coef_perc"] = np.sum(coef_df.to_numpy() == 0, axis=1) / coef_df.shape[1]

    ### Calculate the alpha value where the r2 score significantly dropped
    logging.error("Finding optimal hyperparams...")
    if args.optim_alpha_mode == "kneedle":
        optim_alpha, optim_l1_ratio = optim_alpha_kneedle(args, perf_df)

    ### Save stuff
    perf_df.to_csv(f"{args.save_path}/perf_df.tsv", sep="\t", index=False)
    coef_df.to_csv(f"{args.save_path}/coef_df.tsv", sep="\t", index=False)
    
    return perf_df, coef_df, optim_alpha, optim_l1_ratio


def hyperparam_search_kfold(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):

    alphas, l1_ratios, scores, coeffs, models, conds, pearson_rs, mses, train_inds = [], [], [], [], [], [], [], [], []
    alphas_l = np.logspace(-3, 0, args.num_alphas)
    l1_ratio = 1.0

    # Because here we run kfolds, it's important to do data normalization individually for each training fold
    X_df = atlas_smal_merged.copy()
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[[col, f"log{col}", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    tmp["z_score"] = (2*(tmp[col] > 1) - 1) * tmp["P_value"].apply(lambda x: stats.norm.isf(x / 2))
    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    if args.gene_weight_minmax: weight_col = "-log10(pval)_minmax"
    else: weight_col = "-log10(pval)"

    kf = KFold(n_splits=args.num_folds, shuffle=True, random_state=42)
    data_splits = list(kf.split(tmp))

    for alpha in alphas_l:

        # For each set of params, run kfolds, then decide whether to keep this set of params
        alphas_sub, l1_ratios_sub, scores_sub, coeffs_sub, models_sub, conds_sub, pearson_rs_sub, mses_sub = [], [], [], [], [], [], [], []
        for i, (train_index, test_index) in enumerate(data_splits):

            # logging.error("Hello!!")
            tmp_train, hr_train = get_split_data(tmp, train_index, args.output_label, args.abs_hr)
            tmp_test, hr_test = get_split_data(tmp, test_index, args.output_label, args.abs_hr)

            X_train = tmp_train.copy().drop(columns=[col, f"log{col}", "z_score", "P_value", "-log10(pval)_minmax", "-log10(pval)"])
            X_test = tmp_test.copy().drop(columns=[col, f"log{col}", "z_score", "P_value", "-log10(pval)_minmax", "-log10(pval)"])

            # Normalize the data
            obj = StandardScaler()
            obj = obj.fit(X_train)
            X_train, X_test = obj.transform(X_train), obj.transform(X_test)
    
            # Adjust the alphas carefully, because with full cell-tissue dataset, some alphas do not reach convergence
            model = ElasticNet(l1_ratio=l1_ratio, alpha=alpha, positive=args.positive, fit_intercept=args.intercept, max_iter=5000)

            # Catch whether model converges
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", ConvergenceWarning)
                model.fit(X_train, hr_train)
                if any(issubclass(warning.category, ConvergenceWarning) for warning in w):
                    logging.error(f">>> {args.disease} with alpha={alpha:.3f} and l1_ratio={l1_ratio:.3f} did not converge")
                    continue
        
            # Score the model
            #pred = model.predict(sub_atl)
            #score = r2_score(hr, pred)
            score = model.score(X_test, hr_test)
            pred = model.predict(X_test)
            r, _ = pearsonr(hr_test, pred)
            r = np.nan_to_num(r)
            mse = mean_squared_error(hr_test, pred)
            
            alphas_sub.append(alpha)
            scores_sub.append(score)
            coeffs_sub.append(model.coef_)
            l1_ratios_sub.append(l1_ratio)
            conds_sub.append(args.disease)
            models_sub.append(model)
            pearson_rs_sub.append(r)
            mses_sub.append(mse)
            train_inds.append(train_index)

        # Decide whether to keep this set of params
        alphas.extend(alphas_sub)
        scores.extend(scores_sub)
        coeffs.extend(coeffs_sub)
        l1_ratios.extend(l1_ratios_sub)
        conds.extend(conds_sub)
        models.extend(models_sub)
        pearson_rs.extend(pearson_rs_sub)
        mses.extend(mses_sub)
                
    perf_df = pd.DataFrame({"disease": conds, "alpha": alphas, "l1_ratio": l1_ratios, "pearson_r": pearson_rs, "mse": mses, "r2": scores})
    coef_np = np.array(coeffs)
    # logging.error(f"perf_df.shape: {perf_df.shape}")
    # logging.error(f"coef_np.shape: {coef_np.shape}")
    # logging.error(f"sub_atl.shape: {sub_atl.shape}")
    # logging.error(f"tmp.shape: {tmp.shape}")
    coef_df = pd.DataFrame(coef_np, columns=sub_atl.columns)
    coef_df["disease"] = conds

    # Calculate the average performance across folds
    perf_df["model_run"] = perf_df.apply(lambda x: f"{x['alpha']:.5f}-{x['l1_ratio']:.5f}", axis=1)
    grouped_df = perf_df.groupby("model_run")[["r2", "pearson_r", "mse"]].mean().reset_index().rename(columns={"r2": "r2_mean", "pearson_r": "pearson_r_mean", "mse": "mse_mean"})
    perf_df = perf_df.merge(grouped_df, on="model_run", how="left")

    # Found the set of params with the highest performance in pearson's R and R^2 respectively
    if args.kfold_param == "r2": top_ind = perf_df["r2_mean"].idxmax()
    elif args.kfold_param == "pearson_r": top_ind = perf_df["pearson_r_mean"].idxmax()
    elif args.kfold_param == "mse": top_ind = perf_df["mse_mean"].idxmin()
    top_model = perf_df.loc[top_ind, "model_run"]
    top_alpha, top_l1 = float(top_model.split("-")[0]), float(top_model.split("-")[1])
    top_r2, top_pearson, top_mse = perf_df.loc[top_ind, "r2_mean"], perf_df.loc[top_ind, "pearson_r_mean"], perf_df.loc[top_ind, "mse_mean"]

    with open(f"{args.save_path}/best_model_by_{args.kfold_param}.txt", "w") as f:
        f.write(f"Best model by {args.kfold_param} has alpha={top_alpha:.3f}, l1_ratio={top_l1:.3f} \
                with mean kfold r2={top_r2:.3f}, pearson_r={top_pearson:.3f}, mse={top_mse:.3f} \n")
    
    return perf_df, coef_df, top_alpha, top_l1


def _sampling(args, dataset_size: int):
    # Sampling subsampling: wo replacement, bootstrapping: with replacement
    if args.samp_method == "subsampling":
        inds = [np.random.choice(np.arange(0, dataset_size), size=args.subset_size, replace=False) for _ in range(args.num_iter)]
    elif args.samp_method == "bootstrapping":
        inds = [np.random.choice(np.arange(0, dataset_size), size=args.subset_size, replace=True) for _ in range(args.num_iter)]
    elif args.samp_method == "subsampling_diff_size":
        inds = [
            np.random.choice(np.arange(0, dataset_size), size=int(ss_perc*dataset_size), replace=False) 
            for ss_perc in np.random.uniform(args.lower_subset_size, args.upper_subset_size, size=args.num_iter)
        ]
    return inds


def _stability_analysis_one_alpha(
    args, obj, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame, 
    l1_ratio: float, alpha: float
):
    # Stability analysis
    coeffs_sub = []
    dataset_size = prot_spec_final.shape[0]
    inds = _sampling(args, dataset_size)

    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"

    for ind in inds:

        atlas_smal_subset = atlas_smal_merged.iloc[ind, :]
        sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
        sub_atl = pd.DataFrame(sub_atl, columns=atlas_smal_subset.columns, index=atlas_smal_subset.index)
        tmp = sub_atl.merge(prot_spec_final[[col, f"log{col}", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)

        tmp, hr, weight_col = prep_data(args, tmp, col)
        sub_atl = sub_atl.loc[tmp.index, :]

        # Adjust the alphas carefully, because with full cell-tissue dataset, some alphas do not reach convergence
        model = ElasticNet(l1_ratio=l1_ratio, alpha=alpha, positive=args.positive, fit_intercept=args.intercept, max_iter=5000)

        # Catch whether model converges
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", ConvergenceWarning)
            if args.gene_weight: model.fit(sub_atl, hr, sample_weight=tmp[weight_col].tolist())
            else: model.fit(sub_atl, hr)
            if any(issubclass(warning.category, ConvergenceWarning) for warning in w):
                logging.error("Model not converges")
                continue

        # Score the model
        #pred = model.predict(sub_atl)
        #score = r2_score(hr, pred)
        if args.gene_weight: score = model.score(sub_atl, hr, sample_weight=tmp[weight_col].tolist())
        else: score = model.score(sub_atl, hr)
        pred = model.predict(sub_atl)
        r, _ = pearsonr(hr, pred)
        r = np.nan_to_num(r)

        coeffs_sub.append(model.coef_)

    # Post analysis
    arr = np.sum(np.array(coeffs_sub) != 0, axis=0) / args.num_iter
    arr_pos = np.sum(np.array(coeffs_sub) > 0, axis=0) / args.num_iter
    arr_neg = np.sum(np.array(coeffs_sub) < 0, axis=0) / args.num_iter

    nonzero_coef_df = pd.DataFrame({"nonzero_coef_perc": arr, "pos_coef_perc": arr_pos, "neg_coef_perc": arr_neg, "cell_tissue": sub_atl.columns})
    agg_arr = np.mean(np.abs(np.array(coeffs_sub)), axis=0)
    agg_coef_df = pd.DataFrame({"agg_abs_coef_val": agg_arr, "cell_tissue": sub_atl.columns})
    coef_sen_df = nonzero_coef_df.merge(agg_coef_df, on="cell_tissue", how="left").sort_values(by="pos_coef_perc", ascending=False).reset_index(drop=True)
    
    # Save stuff
    coef_sen_df.to_csv(f"{args.save_path}/coef_stability_select.tsv", sep="\t", index=False)
    with open(f"{args.save_path}/coeffs_sub_arr.npy", "wb") as f:
        np.save(f, np.array(coeffs_sub))
    with open(f"{args.save_path}/inds.pkl", "wb") as f:
        pickle.dump(inds, f)


def stability_selection(args, alpha: list, l1_ratio: float, prot_spec_final: pd.DataFrame, atlas_smal_merged: pd.DataFrame):
    """
    Run stability selection
    """
    obj = StandardScaler()

    # Model and other params
    dataset_size = prot_spec_final.shape[0]
    if args.subset_size == -1: args.subset_size = (dataset_size) // 2

    _stability_analysis_one_alpha(args, obj, atlas_smal_merged, prot_spec_final, l1_ratio, alpha[0])


def main(args):

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
    logging.error("Loading prot data...")
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)
    prot_spec_final.to_csv(f"{args.save_path}/prot_spec_final.tsv", sep="\t", index=False)

    # Convert some argument values to bool
    args.abs_hr = args.abs_hr == 1
    args.positive = args.pos_coef == 1
    args.gene_weight = args.gene_weight == 1
    args.gene_weight_minmax = args.gene_weight_minmax == 1
    args.intercept = args.intercept == 1


    # Train
    if args.optim_alpha_mode != "alpha_list":
        logging.error("Running hyperparam search...")
        if args.optim_alpha_mode == "kfold":
            perf_df, coef_df, optim_alpha, optim_l1_ratio = hyperparam_search_kfold(args, atlas_smal, prot_spec_final)
        else:
            perf_df, coef_df, optim_alpha, optim_l1_ratio = hyperparam_search(args, atlas_smal, prot_spec_final)
        #logging.error("Running permutation importance for all features...")
        #permute_importance(args, prot_spec_final, atlas_smal, alpha=optim_alpha, l1_ratio=optim_l1_ratio)
        logging.error("Running stability selection...")
        stability_selection(args, [optim_alpha], optim_l1_ratio, prot_spec_final, atlas_smal)
        
    else:
        stability_selection(args, args.alpha_list, 1.0, prot_spec_final, atlas_smal)    
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas_path", type=str, help="Atlas path")
    parser.add_argument("--atlas_smal_path", type=str, help="Atlas smal path")
    parser.add_argument("--prot_data_path", type=str, help="Prot path")
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--disease", type=str)
    parser.add_argument("--output_label", type=str, default="HR")

    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to set HR to abs value")
    parser.add_argument("--pos_coef", type=int, default=0, help="Whether to only have positive coefficients in the model")
    parser.add_argument("--gene_weight", type=int, default=1, help="Whether to use -log10(pval) for gene weight")
    parser.add_argument("--gene_weight_minmax", type=int, default=0, help="Whether to minmax gene weight")
    parser.add_argument("--intercept", type=int, default=0, help="Whether to have intercept in elasticnet")
    parser.add_argument("--num_alphas", type=int, default=100, help="Number of alphas for elasticnet")

    parser.add_argument("--subset_size", type=int, default=-1, help="Subset size")
    parser.add_argument("--lower_subset_size", type=float, default=0.5, help="Lower subset size")
    parser.add_argument("--upper_subset_size", type=float, default=0.8, help="Lower subset size")
    parser.add_argument("--num_iter", type=int, default=1000, help="Number of iterations for stability analysis")
    parser.add_argument("--samp_method", type=str, default="subsampling", help="Sampling method")

    parser.add_argument("--optim_alpha_mode", type=str, default="kneedle", help="Method to find optim alpha value. If set to 'alpha_list', then use --alpha_list")
    parser.add_argument("--alpha_list", type=float, nargs="+", default=[0.2, 0.5, 0.7, 0.8, 0.9, 0.95], help="List of alpha values, overriding --optim_alpha_mode")
    parser.add_argument("--select_thres", type=float, default=0.2, help="Threshold to determine stable variables")
    parser.add_argument("--kfold_param", type=str, default="pearson_r", help="kfold param to choose based on")
    parser.add_argument("--num_folds", type=int, default=10, help="number of folds")

    # For permutation
    parser.add_argument("--kfold_n", type=int, default=5, help="Number of k for kfolds")
    parser.add_argument("--n_permute_repeat", type=int, default=30, help="Number of n permutations")

    args = parser.parse_args()
    main(args)