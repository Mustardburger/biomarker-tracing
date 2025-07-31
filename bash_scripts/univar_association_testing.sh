#!/bin/bash

# Load the modules
module load R/4.3.0

# Get the names of params
script_dir="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/scripts/univar_association_testing.R"
atlas_path=$1
pc_df_path=$2
save_path=$3
prot_data_path=$4
num_pcs=$5
mean_covar=$6
dep_var=$7
z_transform=$8

# Version of the code when splitting the data by cell-tissue pair
Rscript $script_dir \
--exp_df $atlas_path \
--pc_df $pc_df_path \
--prot_df $prot_data_path \
--save_path $save_path \
--num_pc $num_pcs \
--mean_covar $mean_covar \
--dep_var $dep_var \
--z_transform $z_transform