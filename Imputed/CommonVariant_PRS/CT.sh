#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --array=1-11
#SBATCH --mem-per-cpu=20G

module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Imputed/CommonVariant_PRS/CT.R ${SLURM_ARRAY_TASK_ID} > CT"${SLURM_ARRAY_TASK_ID}".Rout