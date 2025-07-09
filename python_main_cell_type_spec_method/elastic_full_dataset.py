import pandas as pd
import numpy as np
import os, argparse, pickle

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import KFold

warnings.simplefilter("ignore", RuntimeWarning)

sns.set_theme()
GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"


# A function to clean some code
def get_split_data(df, inds, abs_hr=False):
    X_train = df.iloc[inds, :]
    hr = np.log(X_train["HR"])
    if abs_hr: hr = np.abs(hr)
    return X_train, hr


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
    prot_df["HR"] = prot_df["HR[95%CI]"].apply(lambda x: float(x.split(" ")[0]))

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
def plot_top_coeff_ens(args, disease, coef_df, perf_df, save_path, num_coeffs=10):
    
    retained_inds = perf_df[(perf_df["pearson_r"] >= args.pearson_r_thres) & (perf_df["r2"] >= args.r2_score_thres)].index
    if len(retained_inds.tolist()) == 0:
        # The cutoff is too stringent. Then pick the top 60 percentile
        pearson_r_cutoff = np.percentile(perf_df["pearson_r"].tolist(), 60)
        r2_score_cutoff = np.percentile(perf_df["r2"].tolist(), 60)
        retained_inds = perf_df[(perf_df["pearson_r"] >= pearson_r_cutoff) & (perf_df["r2"] >= r2_score_cutoff)].index
    else:
        pearson_r_cutoff, r2_score_cutoff = args.pearson_r_thres, args.r2_score_thres


    dis_coef = coef_df[coef_df["disease"] == disease].drop(columns=["disease"]).copy().iloc[retained_inds, :].copy()
    melt_raw = pd.melt(dis_coef).rename(columns={"variable": "cell_tissue", "value": "coefficient"})
    logging.error(melt_raw)
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
    sns.boxplot(data=top_coefs, x="coefficient", y="cell_tissue", hue="labels", dodge=False)
    plt.axvline(0.0, color='black', linestyle='--')
    plt.title(f"Coeffs for models with Pearson's R >= {pearson_r_cutoff:.2f} and R2 >= {r2_score_cutoff}")
    plt.savefig(f"{save_path}/{disease}_features.png", bbox_inches="tight", dpi=300)


# Scatterplot for any disease prediction
def scatter_pred(args, X_df: pd.DataFrame, prot_spec_final: pd.DataFrame, perf_df: pd.DataFrame, disease: str, models: list, name="logHR"):

    #prot_spec_final = generate_prot_spec(disease)
    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    hr = sub_atl.merge(prot_spec_final[["HR", "gene"]].set_index("gene"), right_index=True, left_index=True)["HR"]
    hr = np.log(hr)
    if args.abs_hr: hr = np.abs(hr)

    # Select the alpha with the highest r2
    ind = perf_df[perf_df["disease"] == disease].sort_values(by="r2", ascending=False).head(1).index.item()
    model = models[ind]
    alpha = perf_df.iloc[ind]["alpha"].item()
    l1_ratio = perf_df.iloc[ind]["l1_ratio"].item()
    r2w = perf_df.loc[ind, "r2"]
    pred = model.predict(sub_atl)
    
    score, _ = pearsonr(hr, pred)
    rho, _ = spearmanr(hr, pred)
    slope, intercept = np.polyfit(hr, pred, 1)
    x_trend = np.linspace(min(hr), max(hr), 100)
    y_trend = slope * x_trend + intercept
    
    
    hr_df = hr.to_frame(name=f"true_{name}")
    hr_df[f"pred_{name}"] = pred
    minval = hr_df.min()
    maxval = hr_df.max()

    plt.figure(figsize=(10, 10))
    sns.jointplot(data=hr_df, x=f"true_{name}", y=f"pred_{name}", height=10, s=10)
    plt.plot(np.linspace(minval, maxval, 100), np.linspace(minval, maxval, 100), "r--")
    plt.plot(x_trend, y_trend, color="green")
    plt.axvline(x=0.0, color="black")
    plt.axhline(y=0.0, color="black")
    plt.suptitle(f"Alpha={alpha:.3f}, l1_ratio={l1_ratio:.2f}, r={score:.3f}, rho={rho:.3f}, r2w={r2w:.3f}")
    plt.tight_layout()
    plt.savefig(f"{args.save_path}/{disease}_scatterplot.png", bbox_inches="tight", dpi=300)

    return ind, pred


def train(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):

    # Some evidence that this might work. Let's build this pipeline for all diseases
    # diseases = prot_df["Outcome"].unique().tolist()
    obj = StandardScaler()
    X = obj.fit_transform(atlas_smal_merged.to_numpy())
    X_df = pd.DataFrame(X, columns=atlas_smal_merged.columns, index=atlas_smal_merged.index)

    alphas, l1_ratios, scores, coeffs, models, conds, pearson_rs = [], [], [], [], [], [], []
    alphas_l = np.logspace(-3, 0, args.num_alphas)
    l1_ratios_l = [.2, .5, .7, .8, .85, .9, .95, .975, .99, 1]
    num_ens = len(l1_ratios_l) * len(alphas_l)

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[["HR", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    hr = tmp["HR"]
    hr = np.log(hr)
    if args.abs_hr: hr = np.abs(hr)

    if args.gene_weight_minmax: weight_col = "-log10(pval)_minmax"
    else: weight_col = "-log10(pval)"

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
            # if (r < args.pearson_r_thres) or score < args.r2_score_thres: continue
        
            alphas.append(alpha)
            scores.append(score)
            coeffs.append(model.coef_)
            l1_ratios.append(l1_ratio)
            conds.append(args.disease)
            models.append(model)
            pearson_rs.append(r)
            
                
    perf_df = pd.DataFrame({"disease": conds, "alpha": alphas, "l1_ratio": l1_ratios, "pearson_r": pearson_rs, "r2": scores})
    coef_np = np.array(coeffs)
    coef_df = pd.DataFrame(coef_np, columns=sub_atl.columns)
    coef_df["disease"] = conds

    # Save stuff
    perf_df.to_csv(f"{args.save_path}/perf_df.tsv", sep="\t", index=False)
    coef_df.to_csv(f"{args.save_path}/coef_df.tsv", sep="\t", index=False)
    
    return X_df, perf_df, coef_df, models


def main(args):

    # Save the arguments
    os.makedirs(args.save_path, exist_ok=True)
    args_dict = vars(args)
    with open(f"{args.save_path}/cmd_args.json", "w") as json_file:
        json.dump(args_dict, json_file, indent=4)

    # Load in the atlas data
    full_atlas = pd.read_csv(args.atlas_path, sep="\t")
    atlas_smal = pd.read_csv(args.atlas_smal_path, sep="\t").set_index("gene")

    # Load in prot data
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)

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
    X_df, perf_df, coef_df, models = train(args, atlas_smal, prot_spec_final)

    # Make some plots
    plot_top_coeff_ens(args, args.disease, coef_df, perf_df, args.save_path)
    scatter_pred(args, X_df, prot_spec_final, perf_df, args.disease, models)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas_path", type=str, help="Atlas path")
    parser.add_argument("--atlas_smal_path", type=str, help="Atlas smal path")
    parser.add_argument("--prot_data_path", type=str, help="Prot path")
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--disease", type=str)

    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to set HR to abs value")
    parser.add_argument("--pos_coef", type=int, default=0, help="Whether to only have positive coefficients in the model")
    parser.add_argument("--gene_weight", type=int, default=1, help="Whether to use -log10(pval) for gene weight")
    parser.add_argument("--gene_weight_minmax", type=int, default=1, help="Whether to minmax gene weight")
    parser.add_argument("--intercept", type=int, default=0, help="Whether to have intercept in elasticnet")
    parser.add_argument("--num_alphas", type=int, default=100, help="Number of alphas for elasticnet")
    parser.add_argument("--num_folds", type=int, default=10, help="Number of k folds")
    parser.add_argument("--pearson_r_thres", type=float, default=0.1, help="Minimum pearson r on val set")
    parser.add_argument("--r2_score_thres", type=float, default=0.1, help="Minimum r2 score on val set")

    args = parser.parse_args()
    main(args)