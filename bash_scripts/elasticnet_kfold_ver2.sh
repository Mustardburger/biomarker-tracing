#!/bin/bash

# Load the modules
module load anaconda3/latest
source activate scanpy-env

# Get the names of params
script_dir="python_main_cell_type_spec_method/elastic_kfold_ver2.py"
atlas_path=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
num_alphas=$6
num_folds=$7
gene_weight=$8
output_label=$9
abs_hr=${10}

# Version of the code when splitting the data by cell-tissue pair
python $script_dir \
--atlas_path $atlas_path \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--num_alphas $num_alphas \
--num_folds $num_folds \
--gene_weight $gene_weight \
--output_label $output_label \
--abs_hr $abs_hr