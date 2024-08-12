#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=30G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/AoU/JointPRS_Sumstats.R > JointPRS_Sumstats.Rout