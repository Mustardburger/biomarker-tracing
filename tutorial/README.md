# Prerequisites to run this tutorial
There are 2 parts of the tutorial:
1. Calculate seismic's cell type specificity matrix from a single-cell Seurat object (R script).
2. Run the pipeline to link gene's specificity score to gene's summary statistics (python script).
For the first part, as long as you have R installed and with the necessary packages, it should run. However, for the second part, it will only run on the LSF cluster. It will not run on other platforms. Sorry for the inconvenience.

The exact R package versions used in the first part is stored in `biomarker-tracing/r_packages.txt`, whereas the exact python package versions used in the second part is stored in `biomarker-tracing/python_packages.yml`. Please install these packages first before beginning to run the tutorial.

The necessary files are all in the `biomarker-tracing/tutorial` folder. If you would like to use your own summary statistics file instead of `biomarker-tracing/tutorial/sample_disease/Alcoholic_liver_disease.csv`, format the file similar to `biomarker-tracing/sample_data/sample_sumstats/AD_meta_sumstats_all.csv`.

# Part 1: seismic's cell type specificity metric
To calculate seismic's cell type specificity from a Seurat single-cell atlas object, refer to the script at `biomarker-tracing/cell_tissue_specificity/hpa_seismic_cell_type_specificity.R`. An example Seurat object is provided at `biomarker-tracing/tutorial/human_protein_atlas_single_cell_small.rds`, which is a tiny version of the full Human Protein Atlas single-cell RNA-seq object. After running the script successfully, you should get a tsv file that looks similar to that at `biomarker-tracing/tutorial/human_protein_atlas_seismic_cell_spec_matrix.tsv`.

# Part 2: run the pipeline
There are 2 smaller steps to run the pipeline:
* Create a yml file to specify the paths and hyperparameters to run the pipeline. To reproduce the results in this folder, use the yml file at `biomarker-tracing/tutorial/univar_multivar_sample.yml`. ***Please note*** that the hyperparams under the `univariate` section of the yml file are set such that the results from univariate regression matches the results from the original seismicGWAS package. The hyperparams for multivariate sections are less controlled, and you can play around with them to see if any of them yields the most biologically meaningful results. However, most of the results so far show that the rankings provided by univariate results are already meaningful enough.
* Run the pipeline. Use the script at `biomarker-tracing/tutorial/run_pipeline_part2.sh`. ***Please note*** that in this script, there is this part `module load anaconda3/latest; conda activate scanpy-env`. Change this part to fit your conda environment settings.

