import pandas as pd
import numpy as np
import os, argparse, pickle

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import KFold

warnings.simplefilter("ignore", RuntimeWarning)

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"


# A function to clean some code
def get_split_data(df, inds, output_label, abs_hr=False):
    X_train = df.iloc[inds, :]
    hr = X_train[output_label]
    if abs_hr: hr = np.abs(hr)
    return X_train, hr


# A function to plot top coefficients
def _plot(inds, coef_df: pd.DataFrame, disease: str, save_path: str, save_name: str, num_coeffs: int, alpha: float, l1_ratio: float):
    
    dis_coef = coef_df[coef_df["disease"] == disease].drop(columns=["disease"]).copy().iloc[inds, :].copy()
    melt_raw = pd.melt(dis_coef).rename(columns={"variable": "cell_tissue", "value": "coefficient"})
    logging.error(melt_raw.shape)
    melt_sort_df = (
        melt_raw
        .groupby("cell_tissue")
        .mean()
        .reset_index()
        .sort_values(by="coefficient", ascending=False)
    )
    top_cel_tis = melt_sort_df.iloc[list(range(num_coeffs)) + list(range(-num_coeffs, 0))]["cell_tissue"].tolist()
    
    top_coefs = melt_raw[melt_raw["cell_tissue"].isin(top_cel_tis)].copy()
    top_coefs['tmp'] = pd.Categorical(top_coefs['cell_tissue'], categories=top_cel_tis, ordered=True)
    top_coefs = top_coefs.sort_values('tmp')
    top_coefs["labels"] = top_coefs["coefficient"].apply(lambda x: "positive" if x > 0 else "negative")
    
    plt.figure(figsize=(9, 6))
    sns.barplot(data=top_coefs, x="coefficient", y="cell_tissue", hue="labels", dodge=False)
    plt.axvline(0.0, color='black', linestyle='--')
    plt.title(f"Coeffs across kfolds of model with alpha={alpha:.3f} and l1_ratio={l1_ratio:.3f}")
    plt.savefig(f"{save_path}/{save_name}.png", bbox_inches="tight", dpi=300)


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

    risk = "HR[95%CI]"
    if risk not in prot_df.columns: risk = "OR[95%CI]"
    risk_sm = risk.split("[")[0]
    prot_df[risk_sm] = prot_df[risk].apply(lambda x: float(x.split(" ")[0]))
    prot_df[f"log{risk_sm}"] = np.log(prot_df[risk_sm])

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


# For each disease, rank and plot the top most positive and top most negative coefficients
def plot_top_coeff_ens(args, disease: str, coef_df: pd.DataFrame, perf_df: pd.DataFrame, 
                       coef_df_full: pd.DataFrame, full_model_df: pd.DataFrame, 
                       save_path: str, num_coeffs=10):

    # Get the top params
    model_runs = full_model_df["model_run"].unique().tolist()

    # Get the coeffs for individual data folds
    for model_run in model_runs:
        alpha, l1_ratio = float(model_run.split("-")[0]), float(model_run.split("-")[1])
        
        # Retain the top inds for kfolds
        inds = perf_df[(perf_df["model_run"] == model_run)].index
        _plot(inds, coef_df, disease, save_path, f"{disease}_{model_run}_kfolds", num_coeffs, alpha, l1_ratio)

        # Retain the top inds for full model
        inds = full_model_df[(full_model_df["model_run"] == model_run)].index
        _plot(inds, coef_df_full, disease, save_path, f"{disease}_{model_run}_full", num_coeffs, alpha, l1_ratio)


def train(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):

    # Some evidence that this might work. Let's build this pipeline for all diseases
    # diseases = prot_df["Outcome"].unique().tolist()

    alphas, l1_ratios, scores, coeffs, models, conds, pearson_rs, mses, train_inds = [], [], [], [], [], [], [], [], []
    alphas_l = np.logspace(-3, 0, args.num_alphas)
    # l1_ratios_l = [.2, .5, .7, .9, .95, .99, 1]
    l1_ratios_l = [0.001, 0.005, 0.01, 0.05, 0.1, 0.15, .2, .7, .9, 1]
    num_ens = len(l1_ratios_l) * len(alphas_l)

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
    # logging.error(data_splits)

    for alpha in alphas_l:
        for l1_ratio in l1_ratios_l:

            logging.error(f"Working on alpha={alpha:.2f} l1_ratio={l1_ratio:.2f}")
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
                    if args.gene_weight: model.fit(X_train, hr_train, sample_weight=tmp_train[weight_col].tolist())
                    else: model.fit(X_train, hr_train)
                    if any(issubclass(warning.category, ConvergenceWarning) for warning in w):
                        logging.error(f">>> {args.disease} with alpha={alpha:.3f} and l1_ratio={l1_ratio:.3f} did not converge")
                        continue
            
                # Score the model
                #pred = model.predict(sub_atl)
                #score = r2_score(hr, pred)
                if args.gene_weight: score = model.score(X_test, hr_test, sample_weight=tmp_test[weight_col].tolist())
                else: score = model.score(X_test, hr_test)
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
    logging.error(f"perf_df.shape: {perf_df.shape}")
    logging.error(f"coef_np.shape: {coef_np.shape}")
    logging.error(f"sub_atl.shape: {sub_atl.shape}")
    logging.error(f"tmp.shape: {tmp.shape}")
    coef_df = pd.DataFrame(coef_np, columns=sub_atl.columns)
    coef_df["disease"] = conds

    # Calculate the average performance across folds
    perf_df["model_run"] = perf_df.apply(lambda x: f"{x['alpha']:.5f}-{x['l1_ratio']:.5f}", axis=1)
    grouped_df = perf_df.groupby("model_run")[["r2", "pearson_r", "mse"]].mean().reset_index().rename(columns={"r2": "r2_mean", "pearson_r": "pearson_r_mean", "mse": "mse_mean"})
    perf_df = perf_df.merge(grouped_df, on="model_run", how="left")

    # Found the set of params with the highest performance in pearson's R and R^2 respectively
    top_r2_ind = perf_df["r2_mean"].idxmax()
    top_pearson_ind = perf_df["pearson_r_mean"].idxmax()
    top_mse_ind = perf_df["mse_mean"].idxmin()
    top_r2_model = perf_df.loc[top_r2_ind, "model_run"]
    top_pearson_model = perf_df.loc[top_pearson_ind, "model_run"]
    top_mse_model = perf_df.loc[top_mse_ind, "model_run"]

    # Process the full dataset for retraining the top models
    hr = tmp[args.output_label]
    if args.abs_hr: hr = np.abs(hr)
    X_full = tmp.copy().drop(columns=[col, f"log{col}", "z_score", "P_value", "-log10(pval)_minmax", "-log10(pval)"])
    obj = StandardScaler()
    X_full_trans = obj.fit_transform(X_full)

    # Retrain the model using these 2 models on the full dataset
    full_model_df, coeffs_full = [], []
    for model_config in [top_r2_model, top_pearson_model, top_mse_model]:
        alpha, l1_ratio = float(model_config.split("-")[0]), float(model_config.split("-")[1])
        model = ElasticNet(l1_ratio=l1_ratio, alpha=alpha, positive=args.positive, fit_intercept=args.intercept, max_iter=5000)

        # Fit the model
        if args.gene_weight: model.fit(X_full_trans, hr, sample_weight=tmp[weight_col].tolist())
        else: model.fit(X_full_trans, hr)

        # Calculate params
        if args.gene_weight: score = model.score(X_full_trans, hr, sample_weight=tmp[weight_col].tolist())
        else: score = model.score(X_full_trans, hr)
        pred = model.predict(X_full_trans)
        r, _ = pearsonr(hr, pred)
        r = np.nan_to_num(r)
        mse = mean_squared_error(hr, pred)

        # Save results
        full_model_df.append([args.disease, alpha, l1_ratio, r, score, mse, model_config])
        coeffs_full.append(model.coef_)

    full_model_df = pd.DataFrame(full_model_df, columns=["disease", "alpha", "l1_ratio", "pearson_r", "r2", "mse", "model_run"])
    coef_np_full = np.array(coeffs_full)
    coef_df_full = pd.DataFrame(coef_np_full, columns=sub_atl.columns)
    coef_df_full["disease"] = args.disease

    # Save stuff
    perf_df.to_csv(f"{args.save_path}/perf_df.tsv", sep="\t", index=False)
    coef_df.to_csv(f"{args.save_path}/coef_df.tsv", sep="\t", index=False)
    full_model_df.to_csv(f"{args.save_path}/best_model_full_data.tsv", sep="\t", index=False)
    coef_df_full.to_csv(f"{args.save_path}/coef_best_model_full_data.tsv", sep="\t", index=False)
    with open(f"{args.save_path}/train_indices.pkl", "wb") as f:
        pickle.dump(train_inds, f)
    with open(f"{args.save_path}/best_model.txt", "w") as f:
        f.write(f"Best model by r2 has alpha={perf_df.loc[top_r2_ind, 'alpha']:.3f}, l1_ratio={perf_df.loc[top_r2_ind, 'l1_ratio']:.3f} \
                with mean kfold r2={perf_df.loc[top_r2_ind, 'r2_mean']:.3f}, pearson_r={perf_df.loc[top_r2_ind, 'pearson_r_mean']:.3f}, mse={perf_df.loc[top_r2_ind, 'mse_mean']:.3f} \n")
        f.write(f"Best model by pearson's R has alpha={perf_df.loc[top_pearson_ind, 'alpha']:.3f}, l1_ratio={perf_df.loc[top_pearson_ind, 'l1_ratio']:.3f}, \
                with mean kfold r2={perf_df.loc[top_pearson_ind, 'r2_mean']:.3f}, pearson_r={perf_df.loc[top_pearson_ind, 'pearson_r_mean']:.3f}, , mse={perf_df.loc[top_pearson_ind, 'mse_mean']:.3f} \n")
        f.write(f"Best model by MSE has alpha={perf_df.loc[top_mse_ind, 'alpha']:.3f}, l1_ratio={perf_df.loc[top_mse_ind, 'l1_ratio']:.3f}, \
                with mean kfold r2={perf_df.loc[top_mse_ind, 'r2_mean']:.3f}, pearson_r={perf_df.loc[top_mse_ind, 'pearson_r_mean']:.3f}, , mse={perf_df.loc[top_mse_ind, 'mse_mean']:.3f} \n")
    
    return perf_df, coef_df, full_model_df, coef_df_full, num_ens


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
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)
    prot_spec_final.to_csv(f"{args.save_path}/prot_spec_final.tsv", sep="\t", index=False)

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
    perf_df, coef_df, full_model_df, coef_df_full, num_ens = train(args, atlas_smal, prot_spec_final)

    # Make some plots
    plot_top_coeff_ens(args, args.disease, coef_df, perf_df, coef_df_full, full_model_df, args.save_path)
    

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
    parser.add_argument("--num_folds", type=int, default=10, help="Number of k folds")

    args = parser.parse_args()
    main(args)