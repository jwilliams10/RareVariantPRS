#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_Analysis/NullModel.R ${SLURM_ARRAY_TASK_ID} > NullModel"${SLURM_ARRAY_TASK_ID}".Rout