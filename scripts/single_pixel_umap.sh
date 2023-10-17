#!/bin/bash
#SBATCH -J brain_spumap
#SBATCH -A alexandr                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 200G                    # memory pool for all cores
#SBATCH -t 0-30:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o /scratch/trose/slurm.%j.out          # STDOUT
#SBATCH -e /scratch/trose/slurm.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tim.rose@embl.de # send-to address

#module load Anaconda3
cd /home/trose/projects/metaspace_evaluation/slurm_scripts/
#bash conda_slurm.sh
conda activate metabolomics2
python -u single_pixel_umap.py
