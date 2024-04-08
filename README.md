# METASPACE knowledge base analysis

Scripts for all analyses of the METASPACE knowledge to reproduce the Figures of the manuscript:
[Rose et al. "METASPACE: A community-populated knowledge base for spatial metabolomics"](https://doi.org/10.1101/539478)

## Running all analyses

To reproduce all Figures/run all case studies, python scripts and notebooks must be executed in the following order: 

1. Modify the directories for saving all data in the `config.py` file.

2. Download all METASPACE datasets by running `scripts/download_all_datasets.py`

3. Generate general statistics figures of the manuscript (knowledge base statistics): `KB_stats.ipynb`

4. Generate case study 1 of the manuscript (dataset similarity analysis): `CaseStudy1_DatasetSimilarities.ipynb`

5. Create sets of context representative datasets: `get_representative_datasets.ipynb`

6. Download all ion images for representative datasets: `scripts/download_images.py`

7. Compute all colocalizations: `scripts/coloc_pval_{SCENARIO}.py`

8. Perform single pixel integration: `scripts/single_pixel_{SCENARIO}.py`

9. Generate case study 2 of the manuscript (single-pixel analysis): `CaseStudy2_SinglePixel.ipynb` 
(Note that the cluster assignment might change due to random initialization and therefore require manual selection of a new cluster)

10. Generate case study 3 of the manuscript (colocalization analysis): `CaseStudy3_colocalization.ipynb`
(Figure 4 requires the [linex2metaspace](https://github.com/tdrose/lipidranking_metaspace) package, which can be installed from the github repository).

11. Generate plots from case study 4 of the manuscript (co-regulation analysis): Run the R-script in the directory `CaseStudy4` (Further information are included as comments in the top of the file)

12. Create the main figure of the manuscript. This figures uses plots from all case studies: `CaseStudy2_SinglePixel.ipynb` 
(Note that some plots are added to the pdf manually using Inkscape and therefore not added in this notebook).

## Additional notes:
* `scripts/*.py` files usually take a longer time to run and can therefore be submitted to a suitable cluster as jobs (`*.sh` slurm scripts are available for each file).

* The file `figure1_stats.ipynb` requires the table `all_dataset_ids-06-09-23.csv` to create the plot in Figure 1B showing the number of datasets uploaded to METASPACE over time.
Since this table contains IDs of private METASPACE datasets, we cannot make it public.
The code can be commented out to plot only uploaded public datasets over time.

* Results for the publication for the case studies 1-3 have been performed on datasets uploaded before 02.02.23. 
