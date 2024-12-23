#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-120
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/GWAS_SummaryStatistics/SummaryStats_Train.R ${SLURM_ARRAY_TASK_ID} > SummaryStats_Train"${SLURM_ARRAY_TASK_ID}".Rout
