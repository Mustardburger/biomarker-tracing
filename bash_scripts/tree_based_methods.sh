#!/bin/bash

# Load the modules
module load anaconda3/latest
source activate scanpy-env

# Get the names of params
script_dir="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/scripts/tree_based_methods.py"
atlas_path=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
output_label=$6

num_trees=$7
min_samples_split=$8
min_samples_leaf=$9
max_samples=${10}

# Version of the code when splitting the data by cell-tissue pair
python $script_dir \
--atlas_path $atlas_path \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--num_trees $num_trees \
--min_samples_split $min_samples_split \
--min_samples_leaf $min_samples_leaf \
--max_samples $max_samples \
--output_label $output_label