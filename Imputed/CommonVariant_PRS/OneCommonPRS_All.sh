#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=8:00:00
#SBATCH --array=1-11
#SBATCH --mem-per-cpu=10G

module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Imputed/CommonVariant_PRS/OneCommonPRS_All.R ${SLURM_ARRAY_TASK_ID} > OneCommonPRS_All"${SLURM_ARRAY_TASK_ID}".Rout