#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-1879
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/RareVariant_Analysis_Train/STAAR_DynamicWindow.R ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout