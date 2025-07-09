#!/bin/bash

#BSUB -J atlas
#BSUB -P acc_DiseaseGeneCell
#BSUB -n 1
#BSUB -W 1:00
#BSUB -R rusage[mem=4000]
#BSUB -M 30000
#BSUB -eo /sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/logs/%J.stderr 
#BSUB -L /bin/bash

# Load the modules
module load anaconda3/latest
source activate scanpy-env

# Get the names of params
script_dir="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/scripts/lasso_stability_analyses.py"
atlas_path=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
num_alphas=100
subset_size=-1
lower_subset_size=0.6
upper_subset_size=0.8
num_iter=2000
samp_method=subsampling
optim_alpha_mode=alpha_list
alpha_list="0.2 0.3 0.4 0.6 0.7 0.9 1.0"
select_thres=0.1

# Version of the code when splitting the data by cell-tissue pair
python $script_dir \
--atlas_path $atlas_path \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--num_alphas $num_alphas \
--subset_size $subset_size \
--num_iter $num_iter \
--samp_method $samp_method \
--optim_alpha_mode $optim_alpha_mode \
--alpha_list $alpha_list \
--select_thres $select_thres \
--lower_subset_size $lower_subset_size \
--upper_subset_size $upper_subset_size