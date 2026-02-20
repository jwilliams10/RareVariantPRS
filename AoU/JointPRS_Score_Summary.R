# =============================================================================
# [RICE-ANNOTATION] AoU/JointPRS_Score_Summary.R
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
  SCORE_File_META <- NULL
  for(i in 1:22){
    SCORE_File_META <- rbind(SCORE_File_META,read.delim(paste0("/data/williamsjacr/AoU_JointPRS/JointPRS_",trait,"_META_pst_eff_a1_b0.5_phiauto_chr",i,".txt"), header=FALSE))
  }
  write.csv(SCORE_File_META,file = paste0("/data/williamsjacr/AoU_JointPRS/JointPRS_META_Score_",trait,".csv"),row.names = FALSE)
}

