import os, subprocess, argparse

# Instruction to run
# proj_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome"
# python $proj_path/scripts/tree_based_methods_auto-script.py --atlas_smal_path $proj_path/results/atlas_highly_var_genes_merged_corr-thres-0.8_graph-merged.tsv --atlas_path /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/analysis/tabula_sapiens/cell_tissue/specificity_metric/tabula_sapiens_pseudobulk_gene_exp_logcounts.tsv --save_path_suffix random_forest_with_atlas_corr-thres-0.8 --output_label HR


BASH_SCRIPT_DIR = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/bash_scripts/tree_based_methods.sh"


def main(in_args):
    if in_args.output_label == "HR": dis_type = "incident"
    else: dis_type = "prevalent"
    popu_type = "all"
    base_path = f"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data/{dis_type}_popu-{popu_type}"
    diseases = os.listdir(base_path)
    if in_args.save_path_suffix != "": in_args.save_path_suffix = f"_{in_args.save_path_suffix}"

    save_path = f"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/{dis_type}_popu-{popu_type}/tree_based_methods{in_args.save_path_suffix}"

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

        args = [in_args.atlas_path, in_args.atlas_smal_path, base_path, save_full_path, dis_name, in_args.output_label]
        command = ["bsub"] + lsf_params + ["bash", in_args.bash_script_dir] + args
        
        if dis_name == "COPD,_hospital_admissions_1,_only_main_diagnosis":
            subprocess.run(command)
        #subprocess.run(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atlas_smal_path", required=True, type=str)
    parser.add_argument("--atlas_path", required=True, type=str)
    parser.add_argument("--save_path_suffix", required=True, type=str, default="")
    parser.add_argument("--bash_script_dir", required=False, default=BASH_SCRIPT_DIR)
    parser.add_argument("--output_label", required=False, type=str, default="HR")

    args = parser.parse_args()
    main(args)