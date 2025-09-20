import numpy as np
import pandas as pd

def gini_coeff(df, g):
    arr = df.loc[g, :].to_numpy()
    if np.all(arr == 0):
        return np.nan

    arr = arr[arr > 0.0]
    
    x = np.sort(arr)
    n = len(x)
    cumx = np.cumsum(x)
    
    gini = (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    return gini


def main():
    """
    Implement the gini-weighted gene expression matrix
    """
    # Replace data_path and save_path with values of your choice
    data_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/human_protein_atlas_single_cell/cell_tissue/specificity_metric/human_protein_atlas_sc_pseudobulk_gene_exp_logcounts.tsv"
    save_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/atlas_data/human_protein_atlas_all_cell_tissues_gini.tsv"
    
    data = pd.read_csv(data_path, sep="\t").set_index("gene")
    gini_vals = np.array([
        gini_coeff(data, g) for g in data.index.tolist()
    ])
    gene_exp_weighted = data.to_numpy() * gini_vals
    
    final_df = pd.DataFrame(gene_exp_weighted, columns=data.columns, index=data.index).reset_index(names="gene")
    final_df.to_csv(save_path, sep="\t")

if __name__ == "__main__":
    main()