import os, subprocess, argparse

# Instruction to run
# python $proj_path/scripts/univar_association_testing_auto-script.py --atlas_path $proj_path/results/atlas_data/atlas_all_cell_tissues.tsv --save_path_suffix dep-var-z_mean-covar-F_z-trans --z_transform 1
# proj_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome"

BASH_SCRIPT_DIR = "bash_scripts/univar_association_testing.sh"
PC_DF_PATH = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/atlas_data/PC_loadings_plasma_proteome_all_cell_tissues.tsv"


def main(in_args):

    base_path = f"{in_args.disease_prot_dir}/{in_args.disease_type}_popu-{in_args.popu_type}"
    diseases = os.listdir(base_path)
    if in_args.save_path_suffix != "": in_args.save_path_suffix = f"_{in_args.save_path_suffix}"
    save_path = f"{in_args.save_path}/{in_args.disease_type}_popu-{in_args.popu_type}/lasso_stability_analyses{in_args.save_path_suffix}"
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

        args = [in_args.atlas_path, in_args.pc_df_path, save_full_path, in_args.prot_df, str(in_args.num_pcs), 
                in_args.mean_covar, in_args.dep_var, str(in_args.z_transform)]
        command = ["bsub"] + lsf_params + ["-oo", f"{log_path}/{disease}.stdout", "-eo", f"{log_path}/{disease}.stderr"] + ["bash", BASH_SCRIPT_DIR] + args
        
        # if dis_name == "COPD,_hospital_admissions_1,_only_main_diagnosis":
        #    subprocess.run(command)
        subprocess.run(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atlas_smal_path", required=True, type=str)
    parser.add_argument("--atlas_path", required=True, type=str)
    parser.add_argument("--save_path", required=True, type=str)
    parser.add_argument("--save_path_suffix", required=True, type=str, default="")
    parser.add_argument("--prot_df", required=True, type=str, default="")

    parser.add_argument("--disease_prot_dir", required=False, default=DISEASE_PROT_DIR)
    parser.add_argument("--bash_script_dir", required=False, default=BASH_SCRIPT_DIR)
    parser.add_argument("--output_label", required=False, type=str, default="HR")

    parser.add_argument("--pc_df_path", required=False, type=str, default=PC_DF_PATH)
    parser.add_argument("--num_pcs", required=False, type=int, default=0)
    parser.add_argument("--mean_covar", required=False, type=str, default="FALSE")
    parser.add_argument("--dep_var", required=False, type=str, default="z", help="Options: z, beta, neg_log_pval")
    parser.add_argument("--z_transform", required=False, type=int, default=0)
    

    args = parser.parse_args()
    main(args)