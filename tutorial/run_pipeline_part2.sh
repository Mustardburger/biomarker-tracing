#!/bin/bash
#BSUB -P acc_DiseaseGeneCell
#BSUB -J atlas
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -eo "tutorial/log.err"
#BSUB -M 2000
#BSUB -W 1:00

module load anaconda3/latest
conda activate scanpy-env

python pipelines/pipeline_univar_to_multivar.py \
  --yml_file tutorial/univar_multivar_sample.yml
