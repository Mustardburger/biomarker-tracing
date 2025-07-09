import pandas as pd
import numpy as np
import os, argparse, kneed, pickle, logging

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
from sklearn.model_selection import KFold
from sklearn.inspection import permutation_importance

warnings.simplefilter("ignore", RuntimeWarning)

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"
OUTPUT_LABEL = "HR"


# A function to clean some code
def get_split_data(df, inds, abs_hr=False):
    X_train = df.iloc[inds, :]
    hr = np.log(X_train[OUTPUT_LABEL])
    if abs_hr: hr = np.abs(hr)
    return X_train, hr


# Using piecewise function to find optimal alpha
def optim_alpha_piecewise(args, perf_df: pd.DataFrame):
    def piecewise_linear(x, x0, y0, k1, k2):
        return np.piecewise(x, [x < x0], 
                            [lambda x: k1 * x + y0 - k1 * x0, 
                            lambda x: k2 * x + y0 - k2 * x0])

    # Assuming you have your data in x and y arrays
    x = perf_df['zero_coef_perc'].values
    y = perf_df['r2'].values

    # Initial guess for parameters
    p0 = [0.67, -0.2, -0.2, -0.8]  # [x0, y0, k1, k2]

    # Fit the piecewise linear model
    params, _ = curve_fit(piecewise_linear, x, y, p0=p0)
    x0, y0, k1, k2 = params

    # Save the plot that finds the knee point
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y)
    plt.plot(np.sort(x), piecewise_linear(np.sort(x), *params), 'r-', linewidth=2)
    plt.axvline(x=x0, color='g', linestyle='--', label=f'Change point at {x0:.3f}')
    plt.xlabel('zero_coef_perc')
    plt.ylabel('r2')
    plt.title(f"Change point: (x,y)=({x0:.2f},{y0:.2f}), Slopes: {k1:.3f} and {k2:.3f}")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{args.save_path}/piecewise_alpha.png", dpi=300, bbox_inches="tight")

    # Calculate the params at the knee points
    perf_df_ranked = perf_df[(perf_df["zero_coef_perc"] >= x0-0.01) & (perf_df["zero_coef_perc"] <= x0+0.01)].sort_values(by="r2", ascending=False)
    optim_alpha = perf_df_ranked.head(1)["alpha"].item()
    optim_l1_ratio = perf_df_ranked.head(1)["l1_ratio"].item()
    return optim_alpha, optim_l1_ratio


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


# A function to load the proteomics data
def load_prot_data(base_path: str, disease: str, atlas: pd.DataFrame):

    # Mapping gene name to gene ID and vice versa
    gene_names_other = pd.read_csv(GENE_ID_SYMBOLS, sep="\t").rename(columns={"gene_ids": "gene", "gene_symbols": "gene_name"})
    gene_names_hgnc = pd.read_csv(GENE_ID_HGNC, sep="\t")
    gene_names_hgnc = gene_names_hgnc.drop(columns=["Status", "Approved_name", "HGNC_ID", "NCBI_gene_ID", "UCSC_gene_ID"]).rename(columns={"Approved_symbol": "gene_name", "Ensembl_gene_ID": "gene"})

    # Concat 2 gene names dataframe
    gene_names = pd.concat([gene_names_other, gene_names_hgnc], axis=0).drop_duplicates()

    # Some manual mappings from proteomic IDs to gene IDs
    mappings = {
        "BAP18": "ENSG00000258315",  # ID for BACC1, a gene decoding for BAP18
        "CERT": "ENSG00000113163",    # CERT actual name is CERT1 which has ID ENSG00000113163
        "GPR15L": "ENSG00000188373",   # The name on the internet is GPR15LG
        "KIR2DL2": "ENSG00000273661",  # This gene does not have an Ensembl ID that matches in the atlas
        "HLA": "ENSG00000204592",
        "MENT": "ENSG00000143443",  # ID for the homolog C6orf56
        "LEG1": "ENSG00000184530",   # ID for the homolog C6orf58
        "LILRA3": "ENSG00000278046", # This gene does not have an Ensembl ID that matches in the atlas
        "NTproBNP": "ENSG00000120937",
        "HLA-DRA": "ENSG00000204287",    # At row 507, where gene_name is also HLA but predictor is HLA-DRA, the gene ID is ENSG00000204287
        "PALM2": "ENSG00000157654",   # ID for PALM2AKAP2, a fusion gene for PALM2-AKAP2
        "SARG": "ENSG00000182795"   # ID for homolog C1orf116
    }

    # Load in prot data
    prot_df = pd.read_csv(f"{os.path.join(base_path, disease)}.csv")
    prot_df = prot_df.rename(columns={"Protein": "gene_name"})
    prot_df[OUTPUT_LABEL] = prot_df[f"{OUTPUT_LABEL}[95%CI]"].astype(str).apply(lambda x: float(x.split(" ")[0]))

    # prot_spec_id contains some genes with duplicate gene ID
    prot_spec_id = pd.merge(left=prot_df, right=gene_names, on="gene_name", how="left")
    dup_genes = prot_spec_id[prot_spec_id.duplicated('gene_name', keep=False)]["gene"].unique().tolist()

    #### Taking care of genes with missing or duplicate gene IDs ####
    # Find proteins with missing mappings
    a = prot_spec_id[prot_spec_id.duplicated('gene', keep=False)].sort_values(by="gene_name")
    pair = [i.split("_") for i in a["gene_name"].tolist() if "_" in i]
    pair = [item for sublist in pair for item in sublist]
    single = [i for i in a["gene_name"].tolist() if "_" not in i]

    # For singles, the dictionary above provides the mappings
    prot_spec_id['gene'] = prot_spec_id.apply(lambda row: mappings.get(row['gene_name'], row['gene']), axis=1)

    # For pairs, after splitting, most the genes can be mapped to gene IDs
    gene_name_pairs = (
        gene_names[gene_names["gene_name"].isin(pair)]
        .drop_duplicates(subset="gene_name")
        .rename(columns={"gene_name": "gene_name_single"})
    )
    prot_spec_id["gene_name_single"] = prot_spec_id["gene_name"].apply(lambda x: x.split("_")[0] if "_" in x else x)
    prot_spec_id = prot_spec_id.merge(gene_name_pairs, on='gene_name_single', how='left', suffixes=('', '_small'))
    prot_spec_id['gene'] = prot_spec_id['gene'].fillna(prot_spec_id['gene_small'])
    prot_spec_id.drop(columns=['gene_small', "gene_name_single"], inplace=True)

    # After all the mappings, NPPB and NTproBNP has the same ENSG. Retain the one with the smaller pval
    larger_pval = prot_spec_id[prot_spec_id["gene_name"].isin(["NPPB", "NTproBNP"])]["P_value"].max()
    dropped_prot = prot_spec_id[(prot_spec_id["gene_name"].isin(["NPPB", "NTproBNP"])) & (prot_spec_id["P_value"] == larger_pval)]["gene_name"].item()
    prot_spec_id = prot_spec_id[prot_spec_id["gene_name"] != dropped_prot].drop_duplicates(subset="gene")

    #### Also take care of genes with duplicate gene IDs ####
    # Some genes with duplicate gene IDs, the rogue IDs will be omitted when merging with the atlas

    # Some genes are not in the atlas. Let's check what they are
    prot_uniq_genes = set(prot_spec_id["gene"].tolist()) - set(atlas.index.tolist())

    # Remove genes not in atlas
    prot_spec_final = prot_spec_id[~prot_spec_id["gene"].isin(list(prot_uniq_genes))].drop_duplicates(subset="gene", keep="first")
    return prot_spec_final


def hyperparam_search(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):

    # Some evidence that this might work. Let's build this pipeline for all diseases
    # diseases = prot_df["Outcome"].unique().tolist()

    alphas, l1_ratios, scores, coeffs, models, conds, pearson_rs = [], [], [], [], [], [], []
    alphas_l = np.logspace(-2, 0, args.num_alphas)
    l1_ratios_l = [1.0]

    obj = StandardScaler()
    X = obj.fit_transform(atlas_smal_merged.to_numpy())
    X_df = pd.DataFrame(X, columns=atlas_smal_merged.columns, index=atlas_smal_merged.index)

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[[OUTPUT_LABEL, "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    hr = tmp[OUTPUT_LABEL]
    hr = np.log(hr)
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

    ### Calculate the alpha value where the r2 score significantly dropped (2 modes: piecewise function and knee points)
    logging.error("Finding optimal hyperparams...")
    if args.optim_alpha_mode == "piecewise":
        optim_alpha, optim_l1_ratio = optim_alpha_piecewise(args, perf_df)
    elif args.optim_alpha_mode == "kneedle":
        optim_alpha, optim_l1_ratio = optim_alpha_kneedle(args, perf_df)

    ### Save stuff
    perf_df.to_csv(f"{args.save_path}/perf_df.tsv", sep="\t", index=False)
    coef_df.to_csv(f"{args.save_path}/coef_df.tsv", sep="\t", index=False)
    
    return perf_df, coef_df, optim_alpha, optim_l1_ratio


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


def permute_importance(args, prot_spec_final: pd.DataFrame, atlas_smal: pd.DataFrame, alpha: float, l1_ratio: float):
    """
    Calculate permutation importance of features as a negative control
    """
    # Transform the data
    X_df = atlas_smal.copy()

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[[OUTPUT_LABEL, "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    hr = tmp[OUTPUT_LABEL]
    hr = np.log(hr)
    if args.abs_hr: hr = np.abs(hr)
    sub_atl = sub_atl.loc[tmp.index, :]
    feature_names = sub_atl.columns.tolist()

    # Split the data k-fold
    kf = KFold(n_splits=args.kfold_n, shuffle=True)
    df_l = []
    for i, (train_index, test_index) in enumerate(kf.split(sub_atl)):

        obj = StandardScaler()
        X_train, y_train = sub_atl.iloc[train_index], hr.iloc[train_index]
        X_test, y_test = sub_atl.iloc[test_index], hr.iloc[test_index]
        obj.fit(X_train)
        X_train, X_test = obj.transform(X_train), obj.transform(X_test) 

        # Train the model on the train folds
        model = ElasticNet(l1_ratio=l1_ratio, alpha=alpha, positive=args.positive, fit_intercept=args.intercept, max_iter=5000)
        if not args.gene_weight: model.fit(X_train, y_train)
        else: model.fit(X_train, y_train, sample_weight=tmp.iloc[train_index]["-log10(pval)"].tolist())
            
        train_score = model.score(X_train, y_train)
        test_score = model.score(X_test, y_test)
        logging.error(f"Train score without sample weights: {train_score:.3f}, test score without sample weights: {test_score:.3f}")
        if args.gene_weight:
            train_score = model.score(X_train, y_train, sample_weight=tmp.iloc[train_index]["-log10(pval)"].tolist())
            test_score = model.score(X_test, y_test, sample_weight=tmp.iloc[test_index]["-log10(pval)"].tolist())
            logging.error(f"Train score with sample weights: {train_score:.3f}, test score with sample weights: {test_score:.3f}")

        # Validate the model
        if not args.gene_weight: r = permutation_importance(model, X_test, y_test, n_repeats=args.n_permute_repeat, random_state=0)
        else: r = permutation_importance(model, X_test, y_test, sample_weight=tmp.iloc[test_index]["-log10(pval)"].tolist(), n_repeats=args.n_permute_repeat, random_state=0)
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
    plt.title(f"Permutscores at l1_ratio={l1_ratio:.2f}, alpha={alpha:.2f}, kfold={args.kfold_n}, repeats={args.n_permute_repeat}")
    plt.savefig(f"{args.save_path}/permute_importance_scores.png", bbox_inches="tight", dpi=300)


def _stability_analysis_one_alpha(
    args, obj, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame, 
    l1_ratio: float, alpha: float
):
    # Stability analysis
    coeffs_sub = []
    dataset_size = prot_spec_final.shape[0]
    inds = _sampling(args, dataset_size)
    for ind in inds:

        atlas_smal_subset = atlas_smal_merged.iloc[ind, :]
        sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
        sub_atl = pd.DataFrame(sub_atl, columns=atlas_smal_subset.columns, index=atlas_smal_subset.index)
        tmp = sub_atl.merge(prot_spec_final[[OUTPUT_LABEL, "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
        tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
        max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
        tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
        tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
        hr = tmp[OUTPUT_LABEL]
        hr = np.log(hr)
        if args.abs_hr: hr = np.abs(hr)
        sub_atl = sub_atl.loc[tmp.index, :]

        if args.gene_weight_minmax: weight_col = "-log10(pval)_minmax"
        else: weight_col = "-log10(pval)"

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


def _stability_analysis_alpha_list(
    args, obj, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame, 
    l1_ratio: float, alpha_list: list
):
    """
    Calculate a set of stable variables from a list of l1
    """
    summary_arr_l = []
    dataset_size = prot_spec_final.shape[0]
    for alpha in alpha_list: 

        logging.error(f"Working on {alpha:.3f}")
        coeffs_sub = []
        inds = _sampling(args, dataset_size)
        for ind in inds:

            atlas_smal_subset = atlas_smal_merged.iloc[ind, :]
            sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
            sub_atl = pd.DataFrame(sub_atl, columns=atlas_smal_subset.columns, index=atlas_smal_subset.index)
            tmp = sub_atl.merge(prot_spec_final[[OUTPUT_LABEL, "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
            tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
            max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
            tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
            tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
            hr = tmp[OUTPUT_LABEL]
            hr = np.log(hr)
            if args.abs_hr: hr = np.abs(hr)

            if args.gene_weight_minmax: weight_col = "-log10(pval)_minmax"
            else: weight_col = "-log10(pval)"

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

        # Calculate the percentage of params
        arr = np.sum(np.array(coeffs_sub) != 0, axis=0) / args.num_iter
        summary_arr_l.append(arr)

    # Save the full results
    final_arr = np.array(summary_arr_l) # final_arr.shape = (num_alpha, features)
    with open(f"{args.save_path}/coeffs_sub_arr.npy", "wb") as f:
        np.save(f, np.array(final_arr))
    with open(f"{args.save_path}/inds.pkl", "wb") as f:
        pickle.dump(inds, f)

    # Identify the stable variables
    max_vals = np.amax(final_arr, axis=0)
    stable_vars = list(np.where(max_vals >= args.select_thres)[0])

    # Save all results
    final_df = pd.DataFrame(final_arr.T, columns=alpha_list, index=sub_atl.columns)
    final_df["stable"] = [1 if i in stable_vars else 0 for i in range(final_df.shape[0])]
    final_df["max_val"] = max_vals
    final_df.sort_values(by="max_val", ascending=False).to_csv(f"{args.save_path}/coef_stability_select.tsv", sep="\t")


def stability_selection(args, alpha: list, l1_ratio: float, prot_spec_final: pd.DataFrame, atlas_smal_merged: pd.DataFrame):
    """
    Run stability selection
    """
    obj = StandardScaler()

    # Model and other params
    dataset_size = prot_spec_final.shape[0]
    if args.subset_size == -1: args.subset_size = (dataset_size) // 2

    if args.optim_alpha_mode != "alpha_list":
        _stability_analysis_one_alpha(args, obj, atlas_smal_merged, prot_spec_final, l1_ratio, alpha[0])
    else:
        _stability_analysis_alpha_list(args, obj, atlas_smal_merged, prot_spec_final, l1_ratio, alpha)


def main(args):

    # Save the arguments
    global OUTPUT_LABEL
    OUTPUT_LABEL = args.output_label
    os.makedirs(args.save_path, exist_ok=True)
    args_dict = vars(args)
    with open(f"{args.save_path}/cmd_args.json", "w") as json_file:
        json.dump(args_dict, json_file, indent=4)

    # Load in the atlas data
    full_atlas = pd.read_csv(args.atlas_path, sep="\t")
    atlas_smal = pd.read_csv(args.atlas_smal_path, sep="\t").set_index("gene")

    # Load in prot data
    logging.error("Loading prot data...")
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas).dropna(subset=OUTPUT_LABEL)

    # Convert some argument values to bool
    args.abs_hr = args.abs_hr == 1
    args.positive = args.pos_coef == 1
    args.gene_weight = args.gene_weight == 1
    args.gene_weight_minmax = args.gene_weight_minmax == 1
    args.intercept = args.intercept == 1

    # Normalize data
    # obj = StandardScaler()
    # X = obj.fit_transform(atlas_smal.to_numpy())
    # X_df = pd.DataFrame(X, columns=atlas_smal.columns, index=atlas_smal.index)

    # Train
    if args.optim_alpha_mode != "alpha_list":
        logging.error("Running hyperparam search...")
        perf_df, coef_df, optim_alpha, optim_l1_ratio = hyperparam_search(args, atlas_smal, prot_spec_final)
        logging.error("Running permutation importance for all features...")
        permute_importance(args, prot_spec_final, atlas_smal, alpha=optim_alpha, l1_ratio=optim_l1_ratio)
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

    # For permutation
    parser.add_argument("--kfold_n", type=int, default=5, help="Number of k for kfolds")
    parser.add_argument("--n_permute_repeat", type=int, default=30, help="Number of n permutations")

    args = parser.parse_args()
    main(args)