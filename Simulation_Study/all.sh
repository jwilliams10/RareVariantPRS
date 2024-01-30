#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --array=1-200
#SBATCH --mem-per-cpu=30G

# module purge
module load R/4.3.0

## step1:GWAS Summary Statistics from Train Data
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/GWAS_SummaryStatistics/SummaryStats_Train.R ${SLURM_ARRAY_TASK_ID} > SummaryStats_Train"${SLURM_ARRAY_TASK_ID}".Rout
## step2: CT, LDPred,LASSOSUM
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/CommonVariant_PRS/CT.R ${SLURM_ARRAY_TASK_ID} > CT"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/CommonVariant_PRS/LDPred_LASSOSum.R ${SLURM_ARRAY_TASK_ID} > LDPred_LASSOSum"${SLURM_ARRAY_TASK_ID}".Rout &
wait
## step3: SL Common
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/CommonVariant_PRS/OneCommonPRS_All.R ${SLURM_ARRAY_TASK_ID} > OneCommonPRS_All"${SLURM_ARRAY_TASK_ID}".Rout
## step4: Null Models for STAAR
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_Analysis/NullModel.R ${SLURM_ARRAY_TASK_ID} > NullModel"${SLURM_ARRAY_TASK_ID}".Rout
## step5: STAAR Analysis/Effect Sizes/ PRS for sliding window, gene centric coding, gene centric noncoding
## Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_Analysis/STAAR_SlidingWindow.R ${SLURM_ARRAY_TASK_ID} > STAAR_SlidingWindow"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_Analysis/STAAR_GeneCentric_NonCoding.R ${SLURM_ARRAY_TASK_ID} > STAAR_GeneCentric_NonCoding"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_Analysis/STAAR_GeneCentric_Coding.R ${SLURM_ARRAY_TASK_ID} > STAAR_GeneCentric_Coding"${SLURM_ARRAY_TASK_ID}".Rout &
wait
## step6: Rare Variant Standalone PRS and Single Best Rare Variant PRS
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_PRS/RareVariant_Standalone_Best_PRS.R ${SLURM_ARRAY_TASK_ID} > RareVariant_Standalone_Best_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/RareVariant_PRS/Single_RareVariant_PRS_All.R ${SLURM_ARRAY_TASK_ID} > Single_RareVariant_PRS_All"${SLURM_ARRAY_TASK_ID}".Rout &
wait
## step7: Common + Rare PRS
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/Common_Plus_Rare_PRS.R ${SLURM_ARRAY_TASK_ID} > Common_Plus_Rare_PRS"${SLURM_ARRAY_TASK_ID}".Rout