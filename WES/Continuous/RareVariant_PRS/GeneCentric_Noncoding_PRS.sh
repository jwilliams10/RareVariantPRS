#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --array=1-660
#SBATCH --mem-per-cpu=100G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/BMI_GeneCentric_Noncoding_PRS.R ${SLURM_ARRAY_TASK_ID} > BMI_GeneCentric_Noncoding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/LDL_GeneCentric_Noncoding_PRS.R ${SLURM_ARRAY_TASK_ID} > LDL_GeneCentric_Noncoding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/HDL_GeneCentric_Noncoding_PRS.R ${SLURM_ARRAY_TASK_ID} > HDL_GeneCentric_Noncoding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/Height_GeneCentric_Noncoding_PRS.R ${SLURM_ARRAY_TASK_ID} > Height_GeneCentric_Noncoding_PRS"${SLURM_ARRAY_TASK_ID}".Rout & 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/TC_GeneCentric_Noncoding_PRS.R ${SLURM_ARRAY_TASK_ID} > TC_GeneCentric_Noncoding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Continuous/RareVariant_PRS/logTG_GeneCentric_Noncoding_PRS.R ${SLURM_ARRAY_TASK_ID} > logTG_GeneCentric_Noncoding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
wait