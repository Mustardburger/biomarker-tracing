import os, subprocess, argparse

# Instruction to run
# proj_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome"

# Run with alpha_list

BASH_SCRIPT_DIR = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/bash_scripts/univar_regression.sh"

def main(in_args):
    dis_type = "incident"
    popu_type = "all"
    base_path = f"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data/{dis_type}_popu-{popu_type}"
    diseases = os.listdir(base_path)
    if in_args.save_path_suffix != "": in_args.save_path_suffix = f"_{in_args.save_path_suffix}"

    save_path = f"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/{dis_type}_popu-{popu_type}/univar_regression{in_args.save_path_suffix}"

    log_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/logs"
    lsf_params = ["-J", "atlas", "-P", "acc_DiseaseGeneCell", "-n", "1", "-W", "0:30", "-R", "rusage[mem=2000]", "-M", "20000", "-L", "/bin/bash"]
    lsf_params = lsf_params + ["-eo", f"{log_path}/%J.stderr"]

    for disease in sorted(diseases):
        if len(disease.split(".")) == 2:
            dis_name = disease.split(".")[0]
        else:
            dis_name = ".".join(disease.split(".")[:-1])

        dis_name = disease.split(".")[0]
        save_full_path = os.path.join(save_path, dis_name)
        os.makedirs(save_path, exist_ok=True)

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
        command = ["bsub"] + lsf_params + ["bash", in_args.bash_script_dir] + args
        
        #if dis_name == "COPD,_hospital_admissions_1,_only_main_diagnosis":
        #    subprocess.run(command)
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