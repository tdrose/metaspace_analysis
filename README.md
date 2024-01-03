# METASPACE evaluation

Scripts for all analyses of the METASPACE knowledge base and to reproduce the Figures of the manuscript:
[Rose et al. "METASPACE: A community-populated knowledge base for spatial metabolomics"](https://doi.org/10.1101/539478)

## Running all analyses

To reproduce all Figures/run all analyses, python scripts and notebooks must be executed in the following order: 

1. Modify the directories for saving all data in the `config.py` file.

2. Download all METASPACE datasets by running `scripts/download_all_datasets.py`

3. Generate Figure 1 of the manuscript (knowledge base statistics): `figure1_stats.ipynb`

4. Generate Figure 2 of the manuscript (dataset similarity analysis): `figure2_DatasetSimilarities.ipynb`

5. Create sets of context representative datasets: `get_representative_datasets.ipynb`

6. Download all ion images for representative datasets: `scripts/download_images.py`

7. Compute all colocalizations: `scripts/coloc_pval_{SCENARIO}.py`

8. Perform single pixel integration: `scripts/single_pixel_{SCENARIO}.py`

9. Generate Figure 3 of the manuscript (single pixel analysis): `figure3_SinglePixel.ipynb` 
(Note that the cluster assignment might change due to random initialization and therefore require manual selection of a new cluster)

10. Generate Figure 4 of the manuscript (colocalization analysis): `figure4_colocalization.ipynb`
(Figure 4 requires the linex2metaspace package, which is a wrapper for the [linex2 package](https://pypi.org/project/linex2/), which will be made public soon)

11. Generate Figure 5 of the manuscript (co-regulation analysis): **TODO**


## Additional notes:
* `scripts/*.py` files usually take a longer time to run and can therefore might be submitted to a suitable cluster as jobs (`*.sh` slurm scripts are available for each file).

* The file `figure1_stats.ipynb` requires the table `all_dataset_ids-06-09-23.csv` to create the the plot in Figure 1B showing the overview of datasets uploaded to METASPACE.
Since this table contains IDs of private METASPACE datasets, we cannot make it public.
The code can be commented out to plot only the public uploaded datasets

* Results for the publication have been performed on datasets uploaded before 02.01.23 
