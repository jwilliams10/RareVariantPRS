#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=120G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_Analysis/BMI_STAAR_Gene_Centric_Coding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > BMI"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_Analysis/HDL_STAAR_Gene_Centric_Coding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > HDL"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_Analysis/Height_STAAR_Gene_Centric_Coding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > Height"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_Analysis/LDL_STAAR_Gene_Centric_Coding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > LDL"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_Analysis/logTG_STAAR_Gene_Centric_Coding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > logTG"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_Analysis/TC_STAAR_Gene_Centric_Coding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > TC"${SLURM_ARRAY_TASK_ID}".Rout &
wait