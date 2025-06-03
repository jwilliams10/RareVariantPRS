#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --array=1-55
#SBATCH --mem-per-cpu=30G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Imputed/RareVariant_PRS/Sensitivity_Analysis_RICE_RV.R ${SLURM_ARRAY_TASK_ID} > Sensitivity_Analysis_RICE_RV"${SLURM_ARRAY_TASK_ID}".Rout