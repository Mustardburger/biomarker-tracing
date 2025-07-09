#!/bin/bash

#BSUB -J atlas
#BSUB -P acc_DiseaseGeneCell
#BSUB -n 1
#BSUB -W 2:00
#BSUB -R rusage[mem=4000]
#BSUB -M 30000
#BSUB -eo /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/logs/%J.stderr 
#BSUB -L /bin/bash

# Load the modules
module load anaconda3/latest
source activate scanpy-env

# Get the names of params
script_dir="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/scripts/elastic_kfold_ver2.py"
atlas_path=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
num_alphas=100
num_folds=10
pearson_r_thres=0.1
r2_score_thres=0.1

# Version of the code when splitting the data by cell-tissue pair
python $script_dir \
--atlas_path $atlas_path \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--num_alphas $num_alphas \
--num_folds $num_folds \
--pearson_r_thres $pearson_r_thres \
--r2_score_thres $r2_score_thres