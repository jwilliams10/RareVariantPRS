#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=6:00:00
#SBATCH --array=1-5
#SBATCH --mem-per-cpu=10G

module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/CommonVariant_PRS/Extract_Beta_All.R ${SLURM_ARRAY_TASK_ID} > Extract_Beta_All"${SLURM_ARRAY_TASK_ID}".Rout