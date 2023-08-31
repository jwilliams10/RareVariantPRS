#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=100G

module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/CommonVariant_PRS/CT.R ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout