#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-5
#SBATCH --mem-per-cpu=60G

module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/CommonVariant_PRS/LDPred_LASSOSum.R ${SLURM_ARRAY_TASK_ID} > LDPred_LASSOSum"${SLURM_ARRAY_TASK_ID}".Rout