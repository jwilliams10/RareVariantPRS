#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-200
#SBATCH --mem-per-cpu=20G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/BMI_SlidingWindow_EffectSizes.R ${SLURM_ARRAY_TASK_ID} > BMI_SlidingWindow_EffectSizes"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/LDL_SlidingWindow_EffectSizes.R ${SLURM_ARRAY_TASK_ID} > LDL_SlidingWindow_EffectSizes"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/HDL_SlidingWindow_EffectSizes.R ${SLURM_ARRAY_TASK_ID} > HDL_SlidingWindow_EffectSizes"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/Height_SlidingWindow_EffectSizes.R ${SLURM_ARRAY_TASK_ID} > Height_SlidingWindow_EffectSizes"${SLURM_ARRAY_TASK_ID}".Rout & 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/TC_SlidingWindow_EffectSizes.R ${SLURM_ARRAY_TASK_ID} > TC_SlidingWindow_EffectSizes"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/logTG_SlidingWindow_EffectSizes.R ${SLURM_ARRAY_TASK_ID} > logTG_SlidingWindow_EffectSizes"${SLURM_ARRAY_TASK_ID}".Rout &
wait