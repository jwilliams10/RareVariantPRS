#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-6
#SBATCH --mem-per-cpu=50G

module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/CommonVariant_PRS/Extract_Beta_All.R ${SLURM_ARRAY_TASK_ID} > Extract_Beta_All"${SLURM_ARRAY_TASK_ID}".Rout