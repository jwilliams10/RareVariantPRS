#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-660
#SBATCH --mem-per-cpu=100G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/RareVariant_PRS/Train_GeneCentric_Coding_PRS.R ${SLURM_ARRAY_TASK_ID} > /spin1/home/linux/williamsjacr/RareVariantPRS/RareVariant_PRS/out/out_Train_GeneCentric_Coding_PRS_"${SLURM_ARRAY_TASK_ID}".Rout