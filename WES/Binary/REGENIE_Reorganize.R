rm(list = ls())
load("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.RData")

for(trait in c("BMI","TC","logTG","LDL","HDL","Height","Breast","Prostate","CAD","T2D","Asthma")){
  tmp <- phenotype_train[!is.na(phenotype_train[,trait]),c(2,1,3:ncol(phenotype_train))]
  write.table(tmp,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/All_Train_",trait,"_REGENIE.txt"),sep = '\t',row.names = FALSE,quote = FALSE)
}