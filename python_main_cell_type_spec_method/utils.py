import numpy as np
import pandas as pd
import subprocess, time, re, os, logging
import scipy.stats as stats

GENE_ID_SYMBOLS = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/datasets/Alzheimer/CSF_proteomics_AD_onset/gene_id_symbol_hgnc.tsv"


# A function to load the proteomics data
def load_prot_data(base_path: str, disease: str, atlas: pd.DataFrame):

    # Mapping gene name to gene ID and vice versa
    gene_names_other = pd.read_csv(GENE_ID_SYMBOLS, sep="\t").rename(columns={"gene_ids": "gene", "gene_symbols": "gene_name"})
    gene_names_hgnc = pd.read_csv(GENE_ID_HGNC, sep="\t")
    gene_names_hgnc = gene_names_hgnc.drop(columns=["Status", "Approved_name", "HGNC_ID", "NCBI_gene_ID", "UCSC_gene_ID"]).rename(columns={"Approved_symbol": "gene_name", "Ensembl_gene_ID": "gene"})

    # Concat 2 gene names dataframe
    gene_names = pd.concat([gene_names_other, gene_names_hgnc], axis=0).drop_duplicates()

    # Some manual mappings from proteomic IDs to gene IDs
    mappings = {
        "BAP18": "ENSG00000258315",  # ID for BACC1, a gene decoding for BAP18
        "CERT": "ENSG00000113163",    # CERT actual name is CERT1 which has ID ENSG00000113163
        "GPR15L": "ENSG00000188373",   # The name on the internet is GPR15LG
        "KIR2DL2": "ENSG00000273661",  # This gene does not have an Ensembl ID that matches in the atlas
        "HLA": "ENSG00000204592",
        "MENT": "ENSG00000143443",  # ID for the homolog C6orf56
        "LEG1": "ENSG00000184530",   # ID for the homolog C6orf58
        "LILRA3": "ENSG00000278046", # This gene does not have an Ensembl ID that matches in the atlas
        "NTproBNP": "ENSG00000120937",
        "HLA-DRA": "ENSG00000204287",    # At row 507, where gene_name is also HLA but predictor is HLA-DRA, the gene ID is ENSG00000204287
        "PALM2": "ENSG00000157654",   # ID for PALM2AKAP2, a fusion gene for PALM2-AKAP2
        "SARG": "ENSG00000182795"   # ID for homolog C1orf116
    }

    # Load in prot data
    prot_df = pd.read_csv(f"{os.path.join(base_path, disease)}.csv")
    prot_df = prot_df.rename(columns={"Protein": "gene_name"})

    risk = "HR[95%CI]"
    if risk not in prot_df.columns: risk = "OR[95%CI]"
    risk_sm = risk.split("[")[0]
    prot_df[risk_sm] = prot_df[risk].apply(lambda x: float(x.split(" ")[0]))
    prot_df[f"log{risk_sm}"] = np.log(prot_df[risk_sm])

    # prot_spec_id contains some genes with duplicate gene ID
    prot_spec_id = pd.merge(left=prot_df, right=gene_names, on="gene_name", how="left")
    dup_genes = prot_spec_id[prot_spec_id.duplicated('gene_name', keep=False)]["gene"].unique().tolist()

    #### Taking care of genes with missing or duplicate gene IDs ####
    # Find proteins with missing mappings
    a = prot_spec_id[prot_spec_id.duplicated('gene', keep=False)].sort_values(by="gene_name")
    pair = [i.split("_") for i in a["gene_name"].tolist() if "_" in i]
    pair = [item for sublist in pair for item in sublist]
    single = [i for i in a["gene_name"].tolist() if "_" not in i]

    # For singles, the dictionary above provides the mappings
    prot_spec_id['gene'] = prot_spec_id.apply(lambda row: mappings.get(row['gene_name'], row['gene']), axis=1)

    # For pairs, after splitting, most the genes can be mapped to gene IDs
    gene_name_pairs = (
        gene_names[gene_names["gene_name"].isin(pair)]
        .drop_duplicates(subset="gene_name")
        .rename(columns={"gene_name": "gene_name_single"})
    )
    prot_spec_id["gene_name_single"] = prot_spec_id["gene_name"].apply(lambda x: x.split("_")[0] if "_" in x else x)
    prot_spec_id = prot_spec_id.merge(gene_name_pairs, on='gene_name_single', how='left', suffixes=('', '_small'))
    prot_spec_id['gene'] = prot_spec_id['gene'].fillna(prot_spec_id['gene_small'])
    prot_spec_id.drop(columns=['gene_small', "gene_name_single"], inplace=True)

    # After all the mappings, NPPB and NTproBNP has the same ENSG. Retain the one with the smaller pval
    larger_pval = prot_spec_id[prot_spec_id["gene_name"].isin(["NPPB", "NTproBNP"])]["P_value"].max()
    dropped_prot = prot_spec_id[(prot_spec_id["gene_name"].isin(["NPPB", "NTproBNP"])) & (prot_spec_id["P_value"] == larger_pval)]["gene_name"].item()
    prot_spec_id = prot_spec_id[prot_spec_id["gene_name"] != dropped_prot].drop_duplicates(subset="gene")

    #### Also take care of genes with duplicate gene IDs ####
    # Some genes with duplicate gene IDs, the rogue IDs will be omitted when merging with the atlas

    # Some genes are not in the atlas. Let's check what they are
    prot_uniq_genes = set(prot_spec_id["gene"].tolist()) - set(atlas.index.tolist())

    # Remove genes not in atlas
    prot_spec_final = prot_spec_id[~prot_spec_id["gene"].isin(list(prot_uniq_genes))].drop_duplicates(subset="gene", keep="first")
    return prot_spec_final


# A function to prep the data for training
def prep_data(args, df: pd.DataFrame, col: str):
    """
    Add some columns that help in training
    """

    df["-log10(pval)"] = -np.log10(df["P_value"])
    df["z_score"] = (2*(df[col] > 1) - 1) * df["P_value"].apply(lambda x: stats.norm.isf(x / 2))
    max_non_inf = df.loc[df["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    df = df.replace([np.inf, -np.inf], max_non_inf)
    df["-log10(pval)_minmax"] = (df["-log10(pval)"] - df["-log10(pval)"].min()) / (df["-log10(pval)"].max() - df["-log10(pval)"].min())
    
    y = df[args.output_label]
    if args.abs_hr: y = np.abs(y)

    if args.gene_weight_minmax: weight_col = "-log10(pval)_minmax"
    else: weight_col = "-log10(pval)"

    return df, y, weight_col


def submit_job_and_wait(bsub_cmd: list, wait_time=10):
    # Submit job and capture job ID
    result = subprocess.run(bsub_cmd, capture_output=True, text=True)
    m = re.search(r"<(\d+)>", result.stdout)
    jobid = m.group(1)
    logging.error(f"Submitted job {jobid}, waiting...")

    # Sleep for about 10s to wait for the job to get registered
    time.sleep(10)

    # Poll until job finishes
    while True:
        bjobs = subprocess.run(["bjobs", jobid], capture_output=True, text=True)
        lines = bjobs.stdout.strip().splitlines()
        fields = lines[1].split()

        # 0=JOBID, 1=USER, 2=JOB_NAME, 3=STAT
        stat = fields[3]
        if stat == "DONE" or stat == "EXIT":
            logging.error(f"stat={stat}, job {jobid} finished")
            break
        else:
            logging.error(f"stat={stat}, job {jobid} still running or pending")
        time.sleep(wait_time)