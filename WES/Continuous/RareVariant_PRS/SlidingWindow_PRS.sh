#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --array=1-660
#SBATCH --mem-per-cpu=75G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/BMI_SlidingWindow_PRS.R ${SLURM_ARRAY_TASK_ID} > BMI_SlidingWindow_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/LDL_SlidingWindow_PRS.R ${SLURM_ARRAY_TASK_ID} > LDL_SlidingWindow_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/HDL_SlidingWindow_PRS.R ${SLURM_ARRAY_TASK_ID} > HDL_SlidingWindow_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/Height_SlidingWindow_PRS.R ${SLURM_ARRAY_TASK_ID} > Height_SlidingWindow_PRS"${SLURM_ARRAY_TASK_ID}".Rout & 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/TC_SlidingWindow_PRS.R ${SLURM_ARRAY_TASK_ID} > TC_SlidingWindow_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/logTG_SlidingWindow_PRS.R ${SLURM_ARRAY_TASK_ID} > logTG_SlidingWindow_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
wait