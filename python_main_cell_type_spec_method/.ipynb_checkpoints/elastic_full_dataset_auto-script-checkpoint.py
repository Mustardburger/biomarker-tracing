import os, subprocess, argparse

BASH_SCRIPT_DIR = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/bash_scripts/elasticnet_full_dataset.sh"


def main(in_args):
    dis_type = "incident"
    popu_type = "all"
    base_path = f"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data/{dis_type}_popu-{popu_type}"
    diseases = os.listdir(base_path)
    if in_args.save_path_suffix != "": in_args.save_path_suffix = f"_{in_args.save_path_suffix}"

    save_path = f"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/{dis_type}_popu-{popu_type}/elastic_full_dataset{in_args.save_path_suffix}"

    log_path = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/logs"
    lsf_params = ["-J", "atlas", "-P", "acc_DiseaseGeneCell", "-n", "1", "-W", "0:30", "-R", "rusage[mem=2000]", "-M", "20000", "-L", "/bin/bash"]
    lsf_params = lsf_params + ["-oo", f"{log_path}/%J.stdout", "-eo", f"{log_path}/%J.stderr"]

    for disease in sorted(diseases):
        if len(disease.split(".")) == 2:
            dis_name = disease.split(".")[0]
        else:
            dis_name = ".".join(disease.split(".")[:-1])

        dis_name = disease.split(".")[0]
        save_full_path = os.path.join(save_path, dis_name)
        os.makedirs(save_path, exist_ok=True)

        args = [in_args.atlas_path, in_args.atlas_smal_path, base_path, save_full_path, dis_name]
        command = ["bsub"] + lsf_params + ["bash", BASH_SCRIPT_DIR] + args
        
        # if dis_name == "COPD,_hospital_admissions_1,_only_main_diagnosis":
        # subprocess.run(command)
        if not os.path.exists(f"{save_path}/{dis_name}"):
            print(dis_name)
            subprocess.run(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atlas_smal_path", required=True, type=str)
    parser.add_argument("--atlas_path", required=True, type=str)
    parser.add_argument("--save_path_suffix", required=True, type=str, default="")
    args = parser.parse_args()
    main(args)