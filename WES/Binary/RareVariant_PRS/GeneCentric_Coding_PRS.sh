#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-990
#SBATCH --mem-per-cpu=40G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_PRS/Asthma_GeneCentric_Coding_PRS.R ${SLURM_ARRAY_TASK_ID} > Asthma_GeneCentric_Coding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_PRS/CAD_GeneCentric_Coding_PRS.R ${SLURM_ARRAY_TASK_ID} > CAD_GeneCentric_Coding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_PRS/T2D_GeneCentric_Coding_PRS.R ${SLURM_ARRAY_TASK_ID} > T2D_GeneCentric_Coding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_PRS/Breast_GeneCentric_Coding_PRS.R ${SLURM_ARRAY_TASK_ID} > Breast_GeneCentric_Coding_PRS"${SLURM_ARRAY_TASK_ID}".Rout & 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_PRS/Prostate_GeneCentric_Coding_PRS.R ${SLURM_ARRAY_TASK_ID} > Prostate_GeneCentric_Coding_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
wait