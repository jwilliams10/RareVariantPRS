#!/bin/bash --login
# =============================================================================
# [RICE-ANNOTATION] Simulation_Study8/all.sh
# Purpose: Driver shell script to run the pipeline in this directory end-to-end (often via SLURM job arrays).
#
# Paper linkage:
#   - Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
#   - Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
#   - Supplementary Data: simulation sample sizes and related summaries.
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --array=1-1000
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.2

## step1:GWAS Summary Statistics from Train Data
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/GWAS_SummaryStatistics/SummaryStats_Train.R ${SLURM_ARRAY_TASK_ID} > SummaryStats_Train"${SLURM_ARRAY_TASK_ID}".Rout
## step2: CT, LDPred,LASSOSUM
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/CommonVariant_PRS/CT.R ${SLURM_ARRAY_TASK_ID} > CT"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/CommonVariant_PRS/LDPred_LASSOSum.R ${SLURM_ARRAY_TASK_ID} > LDPred_LASSOSum"${SLURM_ARRAY_TASK_ID}".Rout &
wait
## step3: SL Common
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/CommonVariant_PRS/OneCommonPRS_All.R ${SLURM_ARRAY_TASK_ID} > OneCommonPRS_All"${SLURM_ARRAY_TASK_ID}".Rout
## step4: Null Models for STAAR
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/RareVariant_Analysis/NullModel.R ${SLURM_ARRAY_TASK_ID} > NullModel"${SLURM_ARRAY_TASK_ID}".Rout
## step5: STAAR Analysis/Effect Sizes/ PRS for sliding window, gene centric coding, gene centric noncoding
## Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/RareVariant_Analysis/STAAR_SlidingWindow.R ${SLURM_ARRAY_TASK_ID} > STAAR_SlidingWindow"${SLURM_ARRAY_TASK_ID}".Rout &
## Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/RareVariant_Analysis/STAAR_GeneCentric_NonCoding.R ${SLURM_ARRAY_TASK_ID} > STAAR_GeneCentric_NonCoding"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/RareVariant_Analysis/STAAR_GeneCentric_Coding.R ${SLURM_ARRAY_TASK_ID} > STAAR_GeneCentric_Coding"${SLURM_ARRAY_TASK_ID}".Rout &
wait
## step6: Rare Variant Standalone PRS and Single Best Rare Variant PRS
## Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/RareVariant_PRS/RareVariant_Standalone_Best_PRS.R ${SLURM_ARRAY_TASK_ID} > RareVariant_Standalone_Best_PRS"${SLURM_ARRAY_TASK_ID}".Rout &
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/RareVariant_PRS/Single_RareVariant_PRS_All.R ${SLURM_ARRAY_TASK_ID} > Single_RareVariant_PRS_All"${SLURM_ARRAY_TASK_ID}".Rout &
wait
## step7: Common + Rare PRS
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study8/Common_Plus_Rare_PRS.R ${SLURM_ARRAY_TASK_ID} > Common_Plus_Rare_PRS"${SLURM_ARRAY_TASK_ID}".Rout