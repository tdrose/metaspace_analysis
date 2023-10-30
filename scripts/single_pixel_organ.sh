#!/bin/bash
#SBATCH -J sporgan
#SBATCH -A alexandr                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 100G                    # memory pool for all cores
#SBATCH -t 0-30:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH --array=1-8  # Adjust the range according to your datasets
#SBATCH -o /scratch/trose/slurm.%j.out          # STDOUT
#SBATCH -e /scratch/trose/slurm.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tim.rose@embl.de # send-to address

module load Anaconda3
conda_path=$(which conda)
conda_dir=$(dirname $conda_path)
conda_parent_dir=$(dirname $conda_dir)
conda_sh_path="$conda_parent_dir/etc/profile.d/conda.sh"
source $conda_sh_path
#source /g/easyqbuild/x86_64/Rocky/8/haswell/software/Anaconda3/2023.03-1/etc/profile.d/conda.sh
conda activate metabolomics2
conda env list

cd /home/trose/projects/metaspace_evaluation/scripts/

datasets=("Brain_normal" "Brain_bbknn" "Kidney_normal" "Kidney_bbknn" "Liver_normal" "Liver_bbknn" "Lung_normal" "Lung_bbknn")

# Get the dataset name for this job
dataset_name="${datasets[$SLURM_ARRAY_TASK_ID - 1]}"

python3 -u single_pixel_organ.py $dataset_name
