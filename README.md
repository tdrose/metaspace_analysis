# METASPACE evaluation

> Tim Daniel Rose, contact: tim.rose@embl.de

Scripts & notebooks for the large-scale analysis of METASPACE data.

## Execution

To perform all analyses, the scripts need to be executed in the following order:

* `config.py` - Set the directories for saving (temporary) data

Afterwards, the follwing download & analysis scripts/notebooks need to be run. 
`scripts/*.py` files usually take a longer time to run and can therefore be submitted to a cluster as jobs (`*.sh` slurm scripts available for each file).

1. `download_METASPACE.ipynb` - All data required for the analysis can be downloaded here (After the set of dataset objects is saved, alternatively the script `scripts/download_all_datasets.py` can be used to download the datasets).
2. `datasets_statistics.ipynb` - General statistics on datasets
3. `embedding_HMDB.ipynb` - UMAP visualization of HMDB datasets using total (using only total intensities instead, ignoring spatial information)
4. `datasets_representative.ipynb` - Workflow to find representative datasets for tissues and ions.
5. `scripts/download_images.py` - Download image files for all representative datasets.
6. `scripts/coloc_pval_poslip.py`, `scripts/coloc_pval_posmet.py`, `scripts/coloc_pval_neglip.py` & `scripts/coloc_pval_negmet.py` - Compute colocalization metrics and p-values.
7. `datasets_coloc_analysis.ipynb` - Coloc analysis across selected datasets and comparison to prior knowledge.
8. `scripts/enrichment_pos_lip.R` - Perform enrichment analysis for select sets of molecules (metabolic modules).
9. `datasets_module_analysis.ipynb` - Metabolic module analysis.
