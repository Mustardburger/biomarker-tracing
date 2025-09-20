#!/bin/bash
#BSUB -J atlas
#BSUB -P acc_DiseaseGeneCell
#BSUB -n 1
#BSUB -W 0:30
#BSUB -M 20000
#BSUB -R "rusage[mem=2000]"
#BSUB -eo "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/pipeline.stderr"
#BSUB -oo "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/pipeline.stdout"

cd "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/biomarker-tracing"
yml_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/pipeline_yml"

module load anaconda3/latest
source activate scanpy-env
python python_main_cell_type_spec_method/pipeline_univar_to_multivar.py \
    --yml_file "$yml_path/univar_multivar.yml"
