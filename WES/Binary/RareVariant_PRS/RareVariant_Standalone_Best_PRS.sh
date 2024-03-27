#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_PRS/RareVariant_Standalone_Best_PRS.R ${SLURM_ARRAY_TASK_ID} > RareVariant_Standalone_Best_PRS"${SLURM_ARRAY_TASK_ID}".Rout