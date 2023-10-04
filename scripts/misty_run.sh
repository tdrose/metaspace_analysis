#!/bin/bash
#SBATCH -J mistypl
#SBATCH -A alexandr                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 4                        # number of cores
#SBATCH --mem 100G                    # memory pool for all cores
#SBATCH -t 0-11:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH --array=1-5  # Adjust the range according to your datasets
#SBATCH --output=/scratch/trose/misty/log/output_%A_%a.txt
#SBATCH --error=/scratch/trose/misty/log/error_%A_%a.txt
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tim.rose@embl.de # send-to address

module load R/4.2.2-foss-2022b
datasets=("/scratch/trose/misty/data/Brain_2022-05-31_10h46m34s.csv" "/scratch/trose/misty/data/Brain_2022-05-31_10h27m17s.csv" "/scratch/trose/misty/data/Brain_2022-05-30_20h44m19s.csv" "/scratch/trose/misty/data/Brain_2021-11-11_11h49m37s.csv" "/scratch/trose/misty/data/Brain_2020-05-19_21h56m17s.csv")

# Get the dataset name for this job
dataset_name="${datasets[$SLURM_ARRAY_TASK_ID - 1]}"

Rscript misty_run.R $dataset_name /scratch/trose/misty/results

