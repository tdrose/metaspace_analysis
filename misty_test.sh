#!/bin/bash
#SBATCH -J download_all_datasets
#SBATCH -A alexandr                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 5                        # number of cores
#SBATCH --mem 50G                    # memory pool for all cores
#SBATCH -t 0-11:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o /scratch/trose/slurm.misty.out          # STDOUT
#SBATCH -e /scratch/trose/slurm.misty.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tim.rose@embl.de # send-to address

module load R/4.2.2-foss-2022b
#bash conda_slurm.sh
#conda activate metabolomics
Rscript misty_test.R
