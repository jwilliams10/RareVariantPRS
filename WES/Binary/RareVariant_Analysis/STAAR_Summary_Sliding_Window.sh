#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=60G

module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_Analysis/STAAR_Summary_Sliding_Window.R > STAAR_Summary_Sliding_Window.Rout