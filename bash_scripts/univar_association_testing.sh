#!/bin/bash

# Load the modules
module load anaconda3/latest
source activate scanpy-env

# Get the names of params
script_dir="python_main_cell_type_spec_method/univar_association_testing.py"
atlas_path=$1
atlas_smal_path=$2
prot_data_path=$3
save_path=$4
disease=$5
abs_hr=$6
output_label=$7
covar_df=$8
covar_gini=$9
ztransform_type=${10}

# Version of the code when splitting the data by cell-tissue pair
python $script_dir \
--atlas_path $atlas_path \
--atlas_smal_path $atlas_smal_path \
--prot_data_path $prot_data_path \
--save_path $save_path \
--disease $disease \
--abs_hr $abs_hr \
--output_label $output_label \
--covar_df $covar_df \
--covar_gini $covar_gini \
--ztransform_type $ztransform_type