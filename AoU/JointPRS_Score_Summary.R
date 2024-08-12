rm(list = ls())
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  SCORE_File_META <- NULL
  for(i in 1:22){
    SCORE_File_META <- rbind(SCORE_File_META,read.delim(paste0("/data/williamsjacr/AoU_JointPRS/JointPRS_",trait,"_META_pst_eff_a1_b0.5_phiauto_chr",i,".txt"), header=FALSE))
  }
  write.csv(SCORE_File_META,file = paste0("/data/williamsjacr/AoU_JointPRS/JointPRS_META_Score_",trait,".csv"),row.names = FALSE)
}

