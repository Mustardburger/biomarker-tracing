#!/bin/bash

# Some params
proj_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome"
run_name="seismic_spec_multivar_loghr"     # For consistency, ensure $run_name here is similar to inputs/save_path_suffix in the yml files
yml_save_name="pipeline_yml_${run_name}"

# Before running this pipeline, make sure that the yml files are populated by running python pipelines/utils.py
# python pipelines/utils.py \
#     --save_path /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc//projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results \
#     --save_name $yml_save_name \
#     --yml_full_path /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/biomarker-tracing/pipeline_yml/univar_multivar_seismic.yml \
#     --disease_path /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data/incident_popu-all
# After that, run the rest of this pipeline


# Some paths
yml_path="${proj_path}/results/$yml_save_name"
yml_files=$(ls $yml_path)
elas_path="${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/elastic_kfold_ver2"
stab_path="${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/stability_analyses"
rf_path="${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/tree_based_methods_random_forest"
univar_path="${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/univar_association_testing"
log_path="${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/logs"
mkdir -p $log_path

# Loop through the yml files
for yml_file in ${yml_files[@]};
do
    disease="${yml_file%%.*}"
    yml_full_path="${yml_path}/${yml_file}"
    if [[ ! -f "${elas_path}/${disease}/best_model.txt" || ! -f "${stab_path}/${disease}/regularization_path.pdf" || ! -f "${rf_path}/${disease}/permute_importance_scores.png" || ! -f "${univar_path}/${disease}/univar_regression_results.tsv" ]]; then
        echo $disease
        job="
        #BSUB -P acc_DiseaseGeneCell
        #BSUB -J atlas
        #BSUB -n 1
        #BSUB -R \"span[hosts=1]\"
        #BSUB -R \"rusage[mem=2000]\"
        #BSUB -W 1:00
        #BSUB -M 20000
        #BSUB -eo ${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/logs/pipeline_${disease}.err
        #BSUB -oo ${proj_path}/pipeline_results/pipeline_univar_multivar_${run_name}/incident_popu-all/logs/pipeline_${disease}.out
        
        cd \"/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/biomarker-tracing\"
        
        module load anaconda3/latest
        source activate scanpy-env
        python pipelines/pipeline_univar_to_multivar.py \
            --yml_file \"$yml_full_path\"
        "
        echo "$job" | bsub
    fi
done