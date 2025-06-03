#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --array=1-6
#SBATCH --mem-per-cpu=40G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/AoU/RICECV_AoU_CrossPlatform.R ${SLURM_ARRAY_TASK_ID} > RICECV_AoU_CrossPlatform"${SLURM_ARRAY_TASK_ID}".Rout