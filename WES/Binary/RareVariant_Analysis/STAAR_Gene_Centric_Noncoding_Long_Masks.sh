#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-8
#SBATCH --mem-per-cpu=120G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_Analysis/Asthma_STAAR_Gene_Centric_Noncoding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > Asthma_STAAR_Gene_Centric_Noncoding_Long_Masks"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_Analysis/CAD_STAAR_Gene_Centric_Noncoding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > CAD_STAAR_Gene_Centric_Noncoding_Long_Masks"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_Analysis/T2D_STAAR_Gene_Centric_Noncoding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > T2D_STAAR_Gene_Centric_Noncoding_Long_Masks"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_Analysis/Breast_STAAR_Gene_Centric_Noncoding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > Breast_STAAR_Gene_Centric_Noncoding_Long_Masks"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/Binary/RareVariant_Analysis/Prostate_STAAR_Gene_Centric_Noncoding_Long_Masks.R ${SLURM_ARRAY_TASK_ID} > Prostate_STAAR_Gene_Centric_Noncoding_Long_Masks"${SLURM_ARRAY_TASK_ID}".Rout &
wait