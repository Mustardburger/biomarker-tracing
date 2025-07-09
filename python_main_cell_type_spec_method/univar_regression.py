import pandas as pd
import numpy as np
import os, argparse

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns


from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, spearmanr
from sklearn.feature_selection import f_regression

sns.set_theme()

warnings.simplefilter("ignore", RuntimeWarning)

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"


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
    prot_spec_id['gene'] = prot_spec_id.apply(lambda row: mappings.get(row['gene_name'], row['gene']), axis=1) # type: ignore

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


def run_univariate_regression(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):
    """
    Run univariate regression by sklearn
    """
    obj = StandardScaler()
    X = obj.fit_transform(atlas_smal_merged.to_numpy())
    X_df = pd.DataFrame(X, columns=atlas_smal_merged.columns, index=atlas_smal_merged.index)

    sub_atl = X_df.loc[prot_spec_final["gene"].tolist(), :]
    tmp = sub_atl.merge(prot_spec_final[["HR", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp["-log10(pval)"] = -np.log10(tmp["P_value"])
    max_non_inf = tmp.loc[tmp["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    tmp = tmp.replace([np.inf, -np.inf], max_non_inf)
    tmp["-log10(pval)_minmax"] = (tmp["-log10(pval)"] - tmp["-log10(pval)"].min()) / (tmp["-log10(pval)"].max() - tmp["-log10(pval)"].min())
    hr = tmp["HR"]
    hr = np.log(hr)
    if args.abs_hr: hr = np.abs(hr)
    sub_atl = sub_atl.loc[tmp.index, :]

    f_statistics, pvals = f_regression(sub_atl, hr)

    # Make it into a DataFrame
    result_df = pd.DataFrame({"cell_tissue": sub_atl.columns, "f_stat": f_statistics, "pval": pvals})
    result_df["log_pval"] = -np.log10(result_df["pval"])

    # Save the data
    result_df = result_df.sort_values(by="log_pval", ascending=False)
    result_df.to_csv(f"{args.save_path}/univar_results.tsv", sep="\t", index=False)

    # Plot
    return result_df


def scatterplot(args, result_df: pd.DataFrame, df_other: pd.DataFrame, col: str, name_of_other: str):
    """
    Generate the scatterplot.
    df_other should have a column denoted 'cell_tissue'
    """
    merged_df = result_df.merge(df_other[["cell_tissue", col]].copy(), on="cell_tissue", how="left")

    # Compute correlations
    r, pval_pearson = pearsonr(merged_df[col], merged_df['log_pval'])
    rho, pval_spearman = spearmanr(merged_df[col], merged_df['log_pval'])

    # Create the plot
    # plt.figure(figsize=(8, 6))
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.regplot(data=merged_df, x=col, y='log_pval', scatter_kws={'s': 50}, line_kws={'color': 'red'}, ax=ax)

    # Annotate statistics
    textstr = '\n'.join((
        f"Pearson r = {r:.2f} (p = {pval_pearson:.3g})",
        f"Spearman ρ = {rho:.2f} (p = {pval_spearman:.3g})"
    ))
    ax.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
            fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    top_n = 15
    top_log_pval = merged_df.nlargest(top_n, 'log_pval')["cell_tissue"].tolist()
    top_other = merged_df.nlargest(top_n, col)["cell_tissue"].tolist()
    common = list(set(top_log_pval + top_other))
    top_df = merged_df[merged_df["cell_tissue"].isin(common)]

    for _, row in top_df.iterrows():
        ax.text(
            row[col],
            row['log_pval'],
            row['cell_tissue'],  # adjust this to the name column you want to display
            fontsize=7.5,
            ha='left'
        )

    # Labels
    ax.set_xlabel(name_of_other)
    ax.set_ylabel("Log pval, univariate regression")
    ax.set_title("Comparison of univariate and multivariate feature selection")
    fig.tight_layout()
    fig.savefig(f"{args.save_path}/{name_of_other}.png")
    fig.close() # type: ignore


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

    # Run the model
    result_df = run_univariate_regression(args, atlas_smal, prot_spec_final)

    # Whether to generate a scatterplot to compare results with some other complicated feature selection methods
    if args.df_others != "null":
        args.df_others = args.df_others.split(":")
        args.cols_of_df_others = args.cols_of_df_others.split(":")
        args.names_of_df_others = args.names_of_df_others.split(":")
        for df_other, col, name in zip(args.df_others, args.cols_of_df_others, args.names_of_df_others):
            df = pd.read_csv(df_other, sep="\t")
            scatterplot(args, result_df, df, col, name)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas_path", type=str, help="Atlas path")
    parser.add_argument("--atlas_smal_path", type=str, help="Atlas smal path")
    parser.add_argument("--prot_data_path", type=str, help="Prot path")
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--disease", type=str)
    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to set HR to abs value")

    parser.add_argument("--df_others", type=str, default="null", help="")
    parser.add_argument("--cols_of_df_others", type=str, default="null", help="")
    parser.add_argument("--names_of_df_others", type=str, default="null", help="")

    args = parser.parse_args()
    main(args)