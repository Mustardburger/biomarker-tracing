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

warnings.simplefilter("ignore", RuntimeWarning)

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"
OUTPUT_LABEL = "HR"


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
    logging.error(f"1.1: {prot_df.shape}")
    logging.error(f"1.1: {prot_df}")

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

    logging.error(f"1.2: {prot_spec_id.shape}")
    logging.error(f"1.2: {prot_spec_id}")

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

    logging.error(f"1.3: {prot_spec_final.shape}")
    logging.error(f"1.3: {prot_spec_final}")
    return prot_spec_final


def permute_importance(args, prot_spec_final: pd.DataFrame, atlas_smal: pd.DataFrame):
    """
    Calculate permutation importance of features as a negative control
    """
    # Transform the data
    X_df = atlas_smal.copy()

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

        obj = StandardScaler()
        X_train, y_train = sub_atl.iloc[train_index], hr.iloc[train_index]
        X_test, y_test = sub_atl.iloc[test_index], hr.iloc[test_index]
        obj.fit(X_train)
        X_train, X_test = obj.transform(X_train), obj.transform(X_test) 

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
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"
    obj = StandardScaler()

    atlas_smal_subset = atlas_smal_merged
    sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
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
                    results.append(((n, mf, md, ms), rf.oob_score_))

    # Best results
    best_params, best_oob = max(results, key=lambda x: x[1])
    args.num_trees = best_params["n_estimators"]
    args.max_features = best_params["max_features"]
    args.max_depth = best_params["max_depth"]
    args.min_samples_leaf = best_params["min_samples_leaf"]

    return args


def random_forests(args, prot_spec_final: pd.DataFrame, atlas_smal_merged: pd.DataFrame):
    """
    Run random forests
    """
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"
    obj = StandardScaler()

    atlas_smal_subset = atlas_smal_merged
    sub_atl = obj.fit_transform(atlas_smal_subset.to_numpy())
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
    atlas_smal = pd.read_csv(args.atlas_smal_path, sep="\t").set_index("gene")

    # Load in prot data
    logging.error("Loading prot data...")
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)
    logging.error(f"1: {prot_spec_final.shape}")
    logging.error(f"{prot_spec_final.head()}")

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

    args = parser.parse_args()
    main(args)