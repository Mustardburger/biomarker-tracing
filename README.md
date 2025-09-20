# Package installation
The packages used to run the code is stored at `packages.yml`. One way to install these packages is `conda env create -f packages.yml`.

# Data preparation
## Atlas data
* The downloaded atlas datasets are found at `/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/atlas_data/data`.
* Explore the Tabula Sapiens Seurat object (count number of unique cell-tissue pair, number of cells in each pair, ...): use `tabula_sapiens_cell_tissue_specificity/tabula_sapiens_explore_cell_tissue_step1.R`.
* Generate normalized pseudobulk gene expression for each cell-tissue pair: use `tabula_sapiens_cell_tissue_specificity/tabula_sapiens_cell_spec_mean_expression_step2.1.R`. Assume the mean expression matrix is saved at `/your/path/data/mean_exp.tsv`.
* Generate [seismic](https://ylaboratory.github.io/seismic/) cell type specificity score for each cell-tissue pair: use `tabula_sapiens_cell_tissue_specificity/tabula_sapiens_cell_spec_seismic_step2.2.R`. Assume the mean expression matrix is saved at `/your/path/data/seismic.tsv`.
* For `step3.1` and `step3.2`, the functionalities are similar to those for Step 2, except now the single cells are grouped into tissues only, not cell-type-tissues.
* _(Optional)_ Raw mean expression data is highly multicollinear, which can affect feature selection performance of LASSO. For univariate selection methods, this step is not needed. To further merge or select cell-tissue pairs with similar gene expression profiles, check out `tabula_sapiens_cell_tissue_specificity/merge_similar_cell_tissue.ipynb`. Assume the merged expression matrix is saved at `/your/path/data/merged_exp.tsv`.

## Proteomic data
* Plasma proteomics data from [this paper](https://proteome-phenome-atlas.com/) is saved at `/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/plasma_proteome/data`.


# Individual analyses
* All relevant code in `python_main_cell_type_spec_method`.
* Run the `auto-script` version of the .py code that accepts relevant input arguments. Other arguments can be found in the relavant bash files in the `bash_scripts` folder.

* 4 versions of code:
    * `elastic_kfold_ver2`: Run ElasticNet using k-fold cross-validation, find the best performing model, and obtain cell-tissue coefficients.
    * `lasso_stability_analyses`: Run a custom-made algorithm to find a suitable lambda value for LASSO, then run LASSO stability analysis across 2000 subsamples of the original dataset, then obtain the most stable cell-tissues.
    * `tree_based_methods`: Run random forests and obtain cell-tissue importance scores.
    * `univar_association_testing`: Run univariate regression for each cell-tissue, similar to `univar_regression`, but with more parameters included (beta, z-scores, FDR adjusted, ...) and options to add covariates such as mean gene expression, principal components, ...

* An example way of running LASSO stability: ```python python_main_cell_type_spec_method/lasso_stability_analyses_auto-script.py --atlas_smal_path /your/path/data/merged_exp.tsv --atlas_path /your/path/data/mean_exp.tsv --save_path_suffix trial```.


# Pipeline analysis
Pipelines are designed to run the aforementioned individual scripts in a more systematic way. One such pipeline design is to run univariate regression first, identify the most significant features by FDR correction, then run multivariate methods using these features and rank them based on coefficient magnitudes / permutation importance scores. This is implemented in `pipelines/pipeline_univar_to_multivar.py`.

The single input argument passed into the pipeline is a yml file containing all information needed to run the pipeline, such as input/output paths, disease types, and hyperparameters of individual methods. An example yml file is at `univar_multivar.yml`.

To run a pipeline, launch a lightweight interactive job, then run:
```
python pipelines/pipeline_univar_to_multivar.py --yml_file "$yml_path/univar_multivar.yml"
```
The pipeline will run through the methods in the correct order.
