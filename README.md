Main steps:

# Attention!
The codebase is currently under construction, with many full paths scattered around files. Before running, please check if there are any full paths pointing to irrelavant paths yet.

# Data preparation
## Atlas data
* Download full Tabula Sapiens single-cell object (Tabula Sapiens v1 data)
* Explore the Tabula Sapiens Seurat object (count number of unique cell-tissue pair, number of cells in each pair, ...): use `tabula_sapiens_cell_tissue_specificity/tabula_sapiens_explore_cell_tissue_step1.R`.
* Generate normalized pseudobulk gene expression for each cell-tissue pair: use `tabula_sapiens_cell_tissue_specificity/tabula_sapiens_cell_spec_mean_expression_step2.R`. Assume the mean expression matrix is saved at `/your/path/data/mean_exp.tsv`.
* Merge cell-tissue pairs with similar gene expression profiles: use `tabula_sapiens_cell_tissue_specificity/merge_similar_cell_tissue.ipynb`. Assume the merged expression matrix is saved at `/your/path/data/merged_exp.tsv`.

## Proteomic data
* Download summary statistics for plasma proteomic data


# Run the analysis
* All relevant code in `python_main_cell_type_spec_method`
* Run the `auto-script` version of the .py code that accepts relevant input arguments. Other arguments can be found in the relavant bash files in the `bash_scripts` folder.
* 5 versions of code:
    * `elastic_full_dataset`: Run 1000 ElasticNets of various hyperparameters on the full dataset, and obtain cell-tissue coefficients.
    * `elastic_kfold_ver2`: Run ElasticNet using k-fold cross-validation, find the best performing model, and obtain cell-tissue coefficients.
    * `lasso_stability_analyses`: Run a custom-made algorithm to find a suitable lambda value for LASSO, then run LASSO stability analysis across 2000 subsamples of the original dataset, then obtain the most stable cell-tissues.
    * `tree_based_methods`: Run random forests and obtain cell-tissue importance scores.
    * `univar_regression`: Run univariate regression for each cell-tissue and rank cell-tissues based on F statistic.
* An example way of running LASSO stability: `python python_main_cell_type_spec_method/lasso_stability_analyses_auto-script.py --atlas_smal_path /your/path/data/merged_exp.tsv --atlas_path /your/path/data/mean_exp.tsv --save_path_suffix trial`.
