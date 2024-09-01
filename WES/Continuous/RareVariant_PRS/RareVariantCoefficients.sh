#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-6
#SBATCH --mem-per-cpu=30G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/RareVariantCoefficients.R ${SLURM_ARRAY_TASK_ID} > RareVariantCoefficients"${SLURM_ARRAY_TASK_ID}".Rout