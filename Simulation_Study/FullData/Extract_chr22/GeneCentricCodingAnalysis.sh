#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G

module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/GeneCentricCodingAnalysis.R > out_genecentriccoding.Rout