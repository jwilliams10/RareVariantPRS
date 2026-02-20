# =============================================================================
# [RICE-ANNOTATION] AoU/JointPRS_Score_Copy.R
# Purpose: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
#
# Paper linkage:
#   - Manuscript: Results -> All of Us Results (Fig. 7) and AoU-trained PRS evaluated on UKB (Fig. 8).
#   - Supplementary Figures: Supp. Fig. 17–21 (AoU association diagnostics, PRS performance, and portability).
#   - Supplementary Data: key AoU cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list = ls())
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  system(paste0("gsutil cp /data/williamsjacr/AoU_JointPRS/JointPRS_META_Score_",trait,".csv gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/Continuous/JointPRS/Scores/"))
} 
