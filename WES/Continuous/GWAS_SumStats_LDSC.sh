#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --array=1-6
#SBATCH --mem-per-cpu=20G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/GWAS_SumStats_LDSC.R ${SLURM_ARRAY_TASK_ID} > GWAS_SumStats_LDSC"${SLURM_ARRAY_TASK_ID}".Rout