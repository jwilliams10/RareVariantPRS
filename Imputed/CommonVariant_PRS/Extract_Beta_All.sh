#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --array=1-11
#SBATCH --mem-per-cpu=60G

module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Imputed/CommonVariant_PRS/Extract_Beta_All.R ${SLURM_ARRAY_TASK_ID} > Extract_Beta_All"${SLURM_ARRAY_TASK_ID}".Rout