rm(list = ls())

library(dplyr)

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
NRI_Data_Continuous <- NULL
for(trait in continuous_traits){
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_BestPRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
  pheno_validation$y_validation <- NA
  pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
  
  pheno_validation <- pheno_validation[!is.na(pheno_validation[,trait]),]
  
  for(risk in c(0.01,0.05,0.1)){
    truth_HighRisk <- which(pheno_validation$y_validation > quantile(pheno_validation$y_validation,1 - risk))
    CV_Only_HighRisk <- which(pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk))
    CV_HighRisk_RV_HighRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_HighRisk_RV_LowRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,risk)))
    
    truth_LowRisk <- which(pheno_validation$y_validation < quantile(pheno_validation$y_validation,risk))
    CV_Only_LowRisk <- which(pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,risk))
    CV_LowRisk_RV_HighRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,risk)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_LowRisk_RV_LowRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,risk)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,risk)))
    
    tmp <- data.frame(trait = trait, risk = paste0(100*risk,"%"),
                      Number_LowRisk = length(truth_LowRisk),
                      CV_Only_LowRisk = paste0(sum(CV_Only_LowRisk %in% truth_LowRisk),"(",100*round(sum(CV_Only_LowRisk %in% truth_LowRisk)/length(truth_LowRisk),2),"%)"),
                      CV_LowRisk_RV_LowRisk = paste0(sum(CV_LowRisk_RV_LowRisk %in% truth_LowRisk),"(",100*round(sum(CV_LowRisk_RV_LowRisk %in% truth_LowRisk)/length(truth_LowRisk),2),"%)"),
                      CV_LowRisk_RV_HighRisk = paste0(sum(CV_LowRisk_RV_HighRisk %in% truth_LowRisk),"(",100*round(sum(CV_LowRisk_RV_HighRisk %in% truth_LowRisk)/length(truth_LowRisk),2),"%)"),
                      NRI_LowRisk = (sum(CV_LowRisk_RV_LowRisk %in% truth_LowRisk) - sum(CV_LowRisk_RV_HighRisk %in% truth_LowRisk))/length(truth_LowRisk),
                      Number_HighRisk = length(truth_HighRisk),
                      CV_Only_HighRisk = paste0(sum(CV_Only_HighRisk %in% truth_HighRisk),"(",100*round(sum(CV_Only_HighRisk %in% truth_HighRisk)/length(truth_HighRisk),2),"%)"),
                      CV_HighRisk_RV_LowRisk = paste0(sum(CV_HighRisk_RV_LowRisk %in% truth_HighRisk),"(",100*round(sum(CV_HighRisk_RV_LowRisk %in% truth_HighRisk)/length(truth_HighRisk),2),"%)"),
                      CV_HighRisk_RV_HighRisk = paste0(sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk),"(",100*round(sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk)/length(truth_HighRisk),2),"%)"),
                      NRI_HighRisk = (sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk) - sum(CV_HighRisk_RV_LowRisk %in% truth_HighRisk))/length(truth_HighRisk))
    NRI_Data_Continuous <- rbind(NRI_Data_Continuous,tmp)
  } 
}
NRI_Data_Continuous$NRI_Total <- NRI_Data_Continuous$NRI_LowRisk + NRI_Data_Continuous$NRI_HighRisk

binary_traits <- c("Asthma","T2D","CAD","Breast","Prostate")
NRI_Data_Binary <- NULL
for(trait in binary_traits){
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_BestPRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  pheno_validation <- pheno_validation[!is.na(pheno_validation[,trait]),]
  
  
  for(risk in c(0.01,0.05,0.1)){
    truth_HighRisk <- which(pheno_validation[,trait] == 1)
    CV_Only_HighRisk <- which(pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk))
    CV_HighRisk_RV_HighRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_HighRisk_RV_LowRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,risk)))
    
    truth_LowRisk <- which(pheno_validation[,trait] == 0)
    CV_Only_LowRisk <- which(pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,1 - risk))
    CV_LowRisk_RV_HighRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,risk)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_LowRisk_RV_LowRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,risk)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,risk)))
    
    tmp <- data.frame(trait = trait, risk = paste0(100*risk,"%"),
                      Number_LowRisk = length(truth_LowRisk),
                      CV_Only_LowRisk = paste0(sum(CV_Only_LowRisk %in% truth_LowRisk),"(",100*round(sum(CV_Only_LowRisk %in% truth_LowRisk)/length(truth_LowRisk),2),"%)"),
                      CV_LowRisk_RV_LowRisk = paste0(sum(CV_LowRisk_RV_LowRisk %in% truth_LowRisk),"(",100*round(sum(CV_LowRisk_RV_LowRisk %in% truth_LowRisk)/length(truth_LowRisk),2),"%)"),
                      CV_LowRisk_RV_HighRisk = paste0(sum(CV_LowRisk_RV_HighRisk %in% truth_LowRisk),"(",100*round(sum(CV_LowRisk_RV_HighRisk %in% truth_LowRisk)/length(truth_LowRisk),2),"%)"),
                      NRI_LowRisk = (sum(CV_LowRisk_RV_LowRisk %in% truth_LowRisk) - sum(CV_LowRisk_RV_HighRisk %in% truth_LowRisk))/length(truth_LowRisk),
                      Number_HighRisk = length(truth_HighRisk),
                      CV_Only_HighRisk = paste0(sum(CV_Only_HighRisk %in% truth_HighRisk),"(",100*round(sum(CV_Only_HighRisk %in% truth_HighRisk)/length(truth_HighRisk),2),"%)"),
                      CV_HighRisk_RV_LowRisk = paste0(sum(CV_HighRisk_RV_LowRisk %in% truth_HighRisk),"(",100*round(sum(CV_HighRisk_RV_LowRisk %in% truth_HighRisk)/length(truth_HighRisk),2),"%)"),
                      CV_HighRisk_RV_HighRisk = paste0(sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk),"(",100*round(sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk)/length(truth_HighRisk),2),"%)"),
                      NRI_HighRisk = (sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk) - sum(CV_HighRisk_RV_LowRisk %in% truth_HighRisk))/length(truth_HighRisk))
    NRI_Data_Binary <- rbind(NRI_Data_Binary,tmp)
  } 
}

NRI_Data_Binary$NRI_Total <- NRI_Data_Binary$NRI_LowRisk + NRI_Data_Binary$NRI_HighRisk
