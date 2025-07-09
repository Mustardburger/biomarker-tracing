# Before Tabula Sapiens:

Before looking at the Tabula Sapiens atlas and studying all diseases in the UK-Biobank paper, I only studied biomarkers of Alzheimer's (AD) from cerebrospinal fluid proteomic papers. Therefore, I curated a single-cell dataset of AD at `/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/piehl_2022_gazestani_2023_integrative/results/data/piehl_nph_rerun_finer/sct_merged_piehl_nph_brain_subtypes.rds`. The following code files worked on that seurat object:

* `seismic.R`: seismic is a tool for modeling correlation between cell type specificity and GWAS MAGMA score (bioRxiv, `https://www.biorxiv.org/content/10.1101/2024.05.04.592534v1`). This code remodels seismic into taking in proteomic fold change instead of MAGMA score.
* `seismic_uk_biobank_ppp_binary-class.ipynb`: Jupyer notebook for implementing binary classifier of significant biomarkers of Alzheimer's from UK-Biobank plasma proteomic paper.
* `seismic_uk_biobank_ppp_binary-class_other_neuro_diseases.ipynb`: Same as above, but for other neurodegenerative diseases.


# Present
Then I look at the UK-Biobank plasma proteomics paper. Now the proteins in the blood can come from anywhere in the body, while in the CSF they are mostly confined by the blood-brain barrier. So that's why I have to curate a new dataset that includes all human cells. Some of the relevant code:

* `human_atlas_merge_longrun.R`: Code for merging the human brain cell dataset and the downsampled Tabula Sapiens dataset.
* `human_brain_atlas_explore.R`: Code for exploratory data analysis and processing of the Human Cortical Cell Atlas. In the end, 3 brain regions commonly studied for AD were chosen: middle temporal gyrus (MTG), angular gyrus (AG), and dorsolateral prefrontal cortex (DFC).
* `tabula_sapiens_all_cells_explore.R`: Code for exploratory data analysis and processing of the Tabula Sapiens atlas.
* `human_atlas_uk_biobank_ppp.ipynb`: Jupyter notebook for exploring relationships between the newly curated human single-cell atlas and UK-Biobank dataset.