#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=40G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/Simulate_Data/Simulate_Data_Prop_CausalRareVariants.R > Simulate_Data_Prop_CausalRareVariants.Rout