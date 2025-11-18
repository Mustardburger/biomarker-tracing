import pandas as pd
import numpy as np
import os, argparse, kneed, pickle, logging

import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib

import warnings, json, logging
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet, Lasso
from sklearn.pipeline import Pipeline
from sklearn.utils import check_random_state
from stability_selection import StabilitySelection, plot_stability_path
from utils import *


warnings.simplefilter("ignore", RuntimeWarning)


def train(args, atlas_smal_merged: pd.DataFrame, prot_spec_final: pd.DataFrame):
    """
    Feature selection by stability selection
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

    # Initiate some variables
    obj = StandardScaler()
    col = "HR"
    if col not in prot_spec_final.columns: col = "OR"

    # Get the data
    if (args.ztransform_type == 1): sub_atl = obj.fit_transform(atlas_smal_merged.to_numpy())
    elif (args.ztransform_type == 2): sub_atl = obj.fit_transform(atlas_smal_merged.to_numpy().T).T
    else: sub_atl = atlas_smal_merged.to_numpy()
    sub_atl = pd.DataFrame(sub_atl, columns=atlas_smal_merged.columns, index=atlas_smal_merged.index)
    tmp = sub_atl.merge(prot_spec_final[[col, f"log{col}", "gene", "P_value"]].set_index("gene"), right_index=True, left_index=True)
    tmp, y, weight_col = prep_data(args, tmp, col)
    X = sub_atl.loc[tmp.index, :]

    if args.abs_hr == 1: y = y.abs()

    # Train the thing
    logging.error(X)
    logging.error(y)
    base_estimator = Pipeline([
        ('model', Lasso(max_iter=5000))
    ])

    lambda_array = np.logspace(-1, 0, 20)
    selector = StabilitySelection(base_estimator=base_estimator, lambda_name='model__alpha',
                              lambda_grid=lambda_array).fit(X, y)

    # Retrieve all stability scores
    scores = selector.stability_scores_ # shape = [n_features, n_alphas]
    lambdas = [f"{i:.4f}" for i in lambda_array]
    logging.error(scores)
    score_df = pd.DataFrame(scores, columns=lambdas)
    score_df["feature_label"] = X.columns.tolist()
    score_df.to_csv(f"{args.save_path}/feature_scores.tsv", sep="\t", index=False)

    # Make a plot
    fig, ax = plot_stability_path(selector, threshold_highlight=args.thres)
    fig.savefig(f"{args.save_path}/regularization_path.pdf", bbox_inches="tight")

    # Print stuff
    inds_selected = selector.get_support(indices=True, threshold=args.thres)
    selected_features = X.columns[inds_selected].tolist()
    res = pd.DataFrame({"selected_features": selected_features})
    res.to_csv(f"{args.save_path}/selected_features_thres_{args.thres}.tsv", sep="\t", index=False)


def main(args):

    # Save the arguments
    os.makedirs(args.save_path, exist_ok=True)
    args_dict = vars(args)
    with open(f"{args.save_path}/cmd_args.json", "w") as json_file:
        json.dump(args_dict, json_file, indent=4)

    args.abs_hr = args.abs_hr == 1

    # Load in the atlas data
    full_atlas = pd.read_csv(args.atlas_path, sep="\t")
    if "gene" in full_atlas.columns: full_atlas = full_atlas.set_index("gene")
    atlas_smal = pd.read_csv(args.atlas_smal_path, sep="\t").set_index("gene")

    # Load in prot data
    logging.error("Loading prot data...")
    prot_spec_final = load_prot_data(args.prot_data_path, args.disease, full_atlas)
    prot_spec_final.to_csv(f"{args.save_path}/prot_spec_final.tsv", sep="\t", index=False)

    # Train
    train(args, atlas_smal, prot_spec_final)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas_path", type=str, help="Atlas path")
    parser.add_argument("--atlas_smal_path", type=str, help="Atlas smal path")
    parser.add_argument("--prot_data_path", type=str, help="Prot path")
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--disease", type=str)
    parser.add_argument("--output_label", type=str, default="HR")

    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to set HR to abs value")
    parser.add_argument("--gene_weight_minmax", type=int, default=0, help="Whether to minmax gene weight")
    parser.add_argument("--thres", type=float, default=0.8, help="Cutoff to choose values")

    parser.add_argument("--ztransform_type", type=int, default=1, help="Whether to z transform on each cell type (1) or each gene (2)")

    args = parser.parse_args()
    main(args)