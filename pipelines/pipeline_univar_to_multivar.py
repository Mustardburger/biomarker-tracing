import argparse, yaml, os, subprocess, logging
import pandas as pd
from utils import submit_job_and_wait

# Run this script in a lightweight job
def univar_top_features(config, base_save_path: str, dis_name: str):
    """
    Select top univariate features
    """
    # Select the features based on some threshold
    univar_path = f"{base_save_path}/{config['inputs']['disease_type']}_popu-{config['inputs']['popu_type']}/univar_association_testing/{dis_name}"
    res_df = pd.read_csv(f"{univar_path}/univar_regression_results.tsv", sep="\t")
    col = config["univar_to_multivar"]["column_to_choose"]
    thres = config["univar_to_multivar"]["thres"]
    selected = res_df[res_df[col] <= thres].copy()["cell_tissue"].tolist()

    # Create a new atlas_smal based on the selected features
    atlas_smal_path = pd.read_csv(config["inputs"]["atlas_smal_path"], sep="\t")
    atlas_df_sub = atlas_smal_path[selected + ["gene"]]
    atlas_df_sub.to_csv(
        f"{univar_path}/atlas_smal_path_sig_cel_tis_filtered.tsv", 
        sep="\t", index=False
    )
    return f"{univar_path}/atlas_smal_path_sig_cel_tis_filtered.tsv"


def main(args):
    """
    The pipeline from univariate regression to multivariate regression
    """
    # Load the YAML file
    logging.error(f"Loading yml file...")
    with open(args.yml_file, "r") as f:
        config = yaml.safe_load(f)

    # Some preprocessing and parsing
    logging.error(f"Preprocessing and parsing paths...")
    base_path = f"{config['inputs']['disease_prot_dir']}/{config['inputs']['disease_type']}_popu-{config['inputs']['popu_type']}"
    diseases = os.listdir(base_path)
    if config['inputs']["save_path_suffix"] != "": save_path_suffix = f"_{config['inputs']['save_path_suffix']}"
    else: save_path_suffix = config['inputs']["save_path_suffix"]
    base_save_path = f"{config['inputs']['save_path']}{save_path_suffix}"
    os.makedirs(base_save_path, exist_ok=True)

    # Loop through the diseases
    for disease in sorted(diseases):

        # Check if disease is in yml file
        if len(disease.split(".")) == 2:
            dis_name = disease.split(".")[0]
        else:
            dis_name = ".".join(disease.split(".")[:-1])
        dis_name = disease.split(".")[0]
        if dis_name not in config['inputs']["disease_name"]:
            continue
        else:
            logging.error(f"Start running on {dis_name}...")

        # Run univariate regression
        logging.error(">>> Start running univariate regression")
        logging.error(f">>> Univariate params: {config['univariate']}")
        sub_args = [
            "--atlas_smal_path", config['inputs']['atlas_smal_path'], 
            "--atlas_path", config['inputs']['atlas_path'], 
            "--save_path", base_save_path,
            "--save_path_suffix", "",
            "--disease_name", dis_name,
            "--output_label", config['univariate']['output_label'],
            "--abs_hr", str(config['univariate']['abs_hr']),
            "--covar_df", config['univariate']['covar_df'],
            "--covar_gini", str(config['univariate']['covar_gini'])
        ]
        command = ["python", config['constants']["univar_script_path"]] + sub_args
        submit_job_and_wait(command, wait_time=10)

        # Extract the atlas_smal from the significant features found by univariate
        logging.error(">>> Extract significant features from univariate regression...")
        new_atlas_smal = univar_top_features(config, base_save_path, dis_name)

        # Run elasticnet
        if (config["elasticnet_kfold"]["run"] == 1):
            logging.error(">>> Start running elasticnet...")
            logging.error(f">>> Elasticnet params: {config['elasticnet_kfold']}")
            sub_args = [
                "--atlas_smal_path", new_atlas_smal, 
                "--atlas_path", config['inputs']['atlas_path'], 
                "--save_path", base_save_path,
                "--save_path_suffix", "",
                "--disease_name", dis_name,
                "--output_label", config['elasticnet_kfold']['output_label'],
                "--abs_hr", str(config['elasticnet_kfold']['abs_hr']),
                "--num_alpha", str(config['elasticnet_kfold']['num_alpha']),
                "--num_folds", str(config['elasticnet_kfold']['num_folds']),
                "--gene_weight", str(config['elasticnet_kfold']['gene_weight']),
            ]
            command = ["python", config['constants']["enet_script_path"]] + sub_args
            # submit_job_and_wait(command, wait_time=30)
            subprocess.run(command)

        # Run Lasso stability selection
        if (config["stability_selection"]["run"] == 1):
            logging.error(">>> Start running stability selection...")
            logging.error(f">>> Stability selection params: {config['stability_selection']}")
            sub_args = [
                "--atlas_smal_path", new_atlas_smal, 
                "--atlas_path", config['inputs']['atlas_path'], 
                "--save_path", base_save_path,
                "--save_path_suffix", "",
                "--disease_name", dis_name,
                "--output_label", config['stability_selection']['output_label'],
                "--abs_hr", str(config['stability_selection']['abs_hr']),
                "--thres", str(config['stability_selection']['thres'])
            ]
            command = ["python", config['constants']["stab_sele_script_path"]] + sub_args
            # submit_job_and_wait(command, wait_time=10)
            subprocess.run(command)

        # Run random forests
        if (config["random_forest"]["run"] == 1):
            logging.error(">>> Start running random_forest...")
            logging.error(f">>> Random_forest params: {config['random_forest']}")
            sub_args = [
                "--atlas_smal_path", new_atlas_smal, 
                "--atlas_path", config['inputs']['atlas_path'], 
                "--save_path", base_save_path,
                "--save_path_suffix", "",
                "--disease_name", dis_name,
                "--output_label", config['random_forest']['output_label'],
                "--abs_hr", str(config['random_forest']['abs_hr']),
                "--param_search", str(config['random_forest']['param_search']),
                "--num_trees", str(config['random_forest']['num_trees']),
                "--min_samples_split", str(config['random_forest']['min_samples_split']),
                "--min_samples_leaf", str(config['random_forest']['min_samples_leaf']),
                "--max_samples", str(config['random_forest']['max_samples']),
                "--kfold_n", str(config['random_forest']['kfold_n']),
                "--n_permute_repeat", str(config['random_forest']['n_permute_repeat'])
            ]
            command = ["python", config['constants']["rf_script_path"]] + sub_args
            # submit_job_and_wait(command, wait_time=30)
            subprocess.run(command)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--yml_file", required=True, type=str, help="Path to yml file")

    args = parser.parse_args()
    main(args)