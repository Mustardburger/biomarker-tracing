import os, subprocess, argparse

BASH_SCRIPT_DIR = "bash_scripts/stability_analyses.sh"
DISEASE_PROT_DIR = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data"


def main(in_args):

    base_path = f"{in_args.disease_prot_dir}/{in_args.disease_type}_popu-{in_args.popu_type}"
    diseases = os.listdir(base_path)
    if in_args.save_path_suffix != "": in_args.save_path_suffix = f"_{in_args.save_path_suffix}"
    save_path = f"{in_args.save_path}/{in_args.disease_type}_popu-{in_args.popu_type}/stability_analyses{in_args.save_path_suffix}"
    lsf_params = ["-J", "atlas", "-P", "acc_DiseaseGeneCell", "-n", "1", "-W", "1:30", "-R", "rusage[mem=2000]", "-M", "20000", "-L", "/bin/bash"]

    for disease in sorted(diseases):
        if len(disease.split(".")) == 2:
            dis_name = disease.split(".")[0]
        else:
            dis_name = ".".join(disease.split(".")[:-1])

        dis_name = disease.split(".")[0]
        save_full_path = os.path.join(save_path, dis_name)
        os.makedirs(save_path, exist_ok=True)
        log_path = save_full_path

        args = [
            in_args.atlas_path, in_args.atlas_smal_path, base_path, save_full_path, dis_name, in_args.output_label,
            str(in_args.thres), str(in_args.abs_hr), str(in_args.ztransform_type)
        ]
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
    parser.add_argument("--save_path", required=True, type=str)
    parser.add_argument("--save_path_suffix", required=True, type=str, default="")

    parser.add_argument("--disease_prot_dir", required=False, default=DISEASE_PROT_DIR)
    parser.add_argument("--bash_script_dir", required=False, default=BASH_SCRIPT_DIR)
    parser.add_argument("--output_label", required=False, type=str, default="HR")

    parser.add_argument("--disease_name", required=False, type=str, default="")    
    parser.add_argument("--disease_type", required=False, type=str, default="incident")
    parser.add_argument("--popu_type", required=False, type=str, default="all")

    parser.add_argument("--samp_method", required=False, type=str, default="subsampling_diff_size")
    parser.add_argument("--optim_alpha_mode", required=False, type=str, default="kneedle")
    parser.add_argument("--thres", type=float, default=0.8, help="Cutoff to choose values")
    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to take absolute value of output")   
    parser.add_argument("--ztransform_type", type=int, default=1, help="Whether to z transform on each cell type (1) or each gene (2)") 

    args = parser.parse_args()
    main(args)