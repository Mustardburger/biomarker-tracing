#!/bin/bash
bsub -P acc_DiseaseGeneCell -q interactive -n 1 -R "span[hosts=1]" -R rusage[mem=2000] -W 2:00 -M 20000 -Is /bin/bash

cd "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/biomarker-tracing"
yml_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome/results/pipeline_yml"

module load anaconda3/latest
source activate scanpy-env
python pipelines/pipeline_univar_to_multivar.py \
    --yml_file "$yml_path/univar_multivar.yml"
