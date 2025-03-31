#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=40G

module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Imputed/LDPred_LASSOSum_MapCorr.R > LDPred_LASSOSum_MapCorr.Rout