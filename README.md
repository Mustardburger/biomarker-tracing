# Package installation and systems specification
The packages used to run the code is stored at `packages.yml`. One way to install these packages is `conda env create -f packages.yml`.

**This codebase is designed specifically for HPCs using LSF. This codebase will not work on local machines or Slurm.**

# Data preparation
## Atlas data
The atlas data should have the following format:
* The file should be a tsv file, with columns being the cell types or tissues, and the columns being the genes
* If the file has N columns, then the N-1 columns are the names of the cell types.
* The last column should have the name `gene` and contains either gene Ensembl IDs or gene symbols for each gene in the atlas data.
A sample atlas data can be found at `sample_data/sample_atlas_data/human_protein_atlas_all_tissues.tsv`.


## Proteomic summary statistics data
All the proteomic summary statistics data should be stored in a same folder. The naming of each summary statistic file should follow the format <disease_name>.csv. Each csv file should only contain the summary statistics of all proteins in one particular disease or condition, and different diseases are separated in different files. This codebase is designed to primarily deal with sumstats data in the format introduced from [this paper](https://proteome-phenome-atlas.com/). Examples of such sumstats data can be found at `sample_data/sample_sumstats/Alcoholic_liver_disease.csv` and `sample_data/sample_sumstats/Nonalcoholic_fatty_liver_disease.csv`.

If you would like to provide your own summary statistics, then it should follow the following format:
* It should contain at least these 4 columns (not in a particular order): `P_value`, `logOR` (or `logHR`), `OR` (or `HR`), and `gene`. Other columns in the file will be ignored by the code.
* The `gene` column can contain either gene Ensembl IDs or gene symbols, but it should be in the same format as the `gene` column in the atlas data.
An example custom sumstats file is at `sample_data/sample_sumstats/AD_meta_sumstats_all.csv`


# Individual analyses
* All relevant code in `python_main_cell_type_spec_method`.
* Run the `auto-script` version of the .py code that accepts relevant input arguments. Other arguments can be found in the relavant bash files in the `bash_scripts` folder.

* 4 versions of code:
    * `elastic_kfold_ver2`: Run ElasticNet using k-fold cross-validation, find the best performing model, and obtain cell-tissue coefficients.
    * `lasso_stability_analyses`: Run a custom-made algorithm to find a suitable lambda value for LASSO, then run LASSO stability analysis across 2000 subsamples of the original dataset, then obtain the most stable cell-tissues.
    * `tree_based_methods`: Run random forests and obtain cell-tissue importance scores.
    * `univar_association_testing`: Run univariate regression for each cell-tissue, similar to `univar_regression`, but with more parameters included (beta, z-scores, FDR adjusted, ...) and options to add covariates such as mean gene expression, principal components, ...


# Pipeline analysis
Pipelines are designed to run the aforementioned individual scripts in a more systematic way. One such pipeline design is to run univariate regression first, identify the most significant features by FDR correction, then run multivariate methods using these features and rank them based on coefficient magnitudes / permutation importance scores. This is implemented in `pipelines/pipeline_univar_to_multivar.py`.

The single input argument passed into the pipeline is a yml file containing all information needed to run the pipeline, such as input/output paths, disease types, and hyperparameters of individual methods. An example yml file is at `univar_multivar_sample.yml`.

To run a pipeline, launch a lightweight interactive job, then run:
```
python pipelines/pipeline_univar_to_multivar.py --yml_file pipeline_yml/univar_multivar_sample.yml
```
