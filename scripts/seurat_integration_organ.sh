#!/bin/bash
#SBATCH -J srto
#SBATCH -A alexandr                # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 100G                    # memory pool for all cores
#SBATCH -t 0-50:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH --array=1-4  # Adjust the range according to your datasets
#SBATCH -o /scratch/trose/slurmSEURAT.%j.out          # STDOUT
#SBATCH -e /scratch/trose/slurmSEURAT.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=tim.rose@embl.de # send-to address
module load R/4.2.2-foss-2022b

module load Anaconda3
conda_path=$(which conda)
conda_dir=$(dirname $conda_path)
conda_parent_dir=$(dirname $conda_dir)
conda_sh_path="$conda_parent_dir/etc/profile.d/conda.sh"
source $conda_sh_path
conda activate metabolomics2
conda env list

# Loading the R module is messing with the python intepreter variable...
conda_info=$(conda info)
env_path=$(echo "$conda_info" | awk '/active env location/ {print $5}')


cd /home/trose/projects/metaspace_evaluation/scripts/

datasets=("single_pixel_adata_Brain_normal.pickle" "single_pixel_adata_Kidney_normal.pickle" "single_pixel_adata_Liver_normal.pickle" "single_pixel_adata_Lung_normal.pickle")

# Get the dataset name for this job
dataset_name="${datasets[$SLURM_ARRAY_TASK_ID - 1]}"

$env_path/bin/python3 -u seurat_integration_organ.py $dataset_name
