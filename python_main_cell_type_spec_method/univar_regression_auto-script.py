import os, subprocess, argparse

BASH_SCRIPT_DIR = "bash_scripts/univar_regression.sh"
DISEASE_PROT_DIR = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data"


def main(in_args):
    
    base_path = f"{in_args.disease_prot_dir}/{in_args.dis_type}_popu-{in_args.popu_type}"
    diseases = os.listdir(base_path)
    if in_args.save_path_suffix != "": in_args.save_path_suffix = f"_{in_args.save_path_suffix}"
    save_path = f"{in_args.save_path}/{in_args.dis_type}_popu-{in_args.popu_type}/elastic_kfold_ver2{in_args.save_path_suffix}"
    lsf_params = ["-J", "atlas", "-P", "acc_DiseaseGeneCell", "-n", "1", "-W", "1:30", "-R", "rusage[mem=2000]", "-M", "20000", "-L", "/bin/bash"]
    lsf_params = lsf_params + ["-eo", f"{log_path}/%J.stderr"]

    for disease in sorted(diseases):
        if len(disease.split(".")) == 2:
            dis_name = disease.split(".")[0]
        else:
            dis_name = ".".join(disease.split(".")[:-1])

        dis_name = disease.split(".")[0]
        save_full_path = os.path.join(save_path, dis_name)
        os.makedirs(save_path, exist_ok=True)
        log_path = save_full_path

        if len(in_args.df_other_paths) > 0:
            full_paths = [f"{p}/{dis_name}/{f}" for p, f in zip(in_args.df_other_paths, in_args.df_other_files)]
            full_paths_str = ":".join(full_paths)
            cols_str = ":".join(in_args.cols)
            names_str = ":".join(in_args.names_of_df_other)
        else:
            full_paths_str = "null"
            cols_str = "null"
            names_str = "null"

        args = [in_args.atlas_path, in_args.atlas_smal_path, base_path, save_full_path, dis_name, full_paths_str, cols_str, names_str]
        command = ["bsub"] + lsf_params + ["-oo", f"{log_path}/{disease}.stdout", "-eo", f"{log_path}/{disease}.stderr"] + ["bash", BASH_SCRIPT_DIR] + args
        
        if dis_name == in_args.disease_name:
            subprocess.run(command)
            break
        elif in_args.disease_name == "":
            subprocess.run(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atlas_smal_path", required=True, type=str)
    parser.add_argument("--atlas_path", required=True, type=str)
    parser.add_argument("--save_path_suffix", required=True, type=str, default="")
    parser.add_argument("--bash_script_dir", required=False, default=BASH_SCRIPT_DIR)

    parser.add_argument("--df_other_paths", type=str, nargs="+", default=[])
    parser.add_argument("--df_other_files", type=str, nargs="+", default=[])
    parser.add_argument("--cols", type=str, nargs="+", default=[])
    parser.add_argument("--names_of_df_other", type=str, nargs="+", default=[])

    args = parser.parse_args()
    main(args)