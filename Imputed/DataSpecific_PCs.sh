#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Imputed/DataSpecific_PCs.R ${SLURM_ARRAY_TASK_ID} > DataSpecific_PCs"${SLURM_ARRAY_TASK_ID}".Rout