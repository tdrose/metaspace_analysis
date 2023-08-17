#!/bin/bash
#SBATCH -J poslip_coloc
#SBATCH -A alexandr                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 10                        # number of cores
#SBATCH --mem 250G                    # memory pool for all cores
#SBATCH -t 0-200:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o /scratch/trose/slurm.%j.out          # STDOUT
#SBATCH -e /scratch/trose/slurm.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tim.rose@embl.de # send-to address

#module load Anaconda3
cd /home/trose/projects/metaspace_evaluation/slurm_scripts/
#bash conda_slurm.sh
#conda activate metabolomics
python -u coloc_pval_poslip_combined.py
