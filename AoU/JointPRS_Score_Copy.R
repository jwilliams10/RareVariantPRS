rm(list = ls())
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  system(paste0("gsutil cp /data/williamsjacr/AoU_JointPRS/JointPRS_META_Score_",trait,".csv gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/Continuous/JointPRS/Scores/"))
} 
