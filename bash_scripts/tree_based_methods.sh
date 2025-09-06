#!/bin/bash

# Load the modules
module load anaconda3/latest
source activate scanpy-env

# Get the names of params
script_dir="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/biomarker-tracing/python_main_cell_type_spec_method/tree_based_methods.py"
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

kfold_n=${11}
n_permute_repeat=${12}
param_search=${13}

# Version of the code when splitting the data by cell-tissue pair
python $script_dir \
--atlas_path $atlas_path \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--param_search $param_search \
--num_trees $num_trees \
--min_samples_split $min_samples_split \
--min_samples_leaf $min_samples_leaf \
--max_samples $max_samples \
--kfold_n $kfold_n \
--n_permute_repeat $n_permute_repeat \
--output_label $output_label