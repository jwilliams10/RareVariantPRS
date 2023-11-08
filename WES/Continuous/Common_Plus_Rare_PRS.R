rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

arrayid <- as.numeric(commandArgs(TRUE)[1])

if(arrayid == 1){
  ## STAARO
  
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
  colnames(pheno_tune) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Tune_All.txt")
  
  rarevariant_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_STAARO_Tune_All.txt")
  
  pheno_tune <- left_join(pheno_tune,common_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:16],"CV_PRS")
  pheno_tune <- left_join(pheno_tune,rarevariant_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:17],"RV_PRS")
  
  PRSs_Tune <- pheno_tune[,c(17,18)]
  
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tune)
  
  PRSs_Tune$Residuals <- model.null$residuals
  
  model.best <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
  colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt")
  
  rarevariant_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_STAARO_Validation_All.txt")
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:16],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:17],"RV_PRS")
  
  PRSs_Validation <- pheno_vad[,c(17,18)]
  
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad)
  
  PRSs_Validation$Residuals <- model.null$residuals
  
  predicted_prs <- predict(model.best,PRSs_Validation)
  
  r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square
  
  effects <- summary(model.best)$coefficients[-1,1]
  
  save(effects,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Coefficients_All_STAARO.RData")
  
  data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "CV_plus_RV_STAARO",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/STAARO_All_Result.RData") 
}else{
  ## Burden
  
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
  colnames(pheno_tune) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Tune_All.txt")
  
  rarevariant_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_Burden_Tune_All.txt")
  
  pheno_tune <- left_join(pheno_tune,common_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:16],"CV_PRS")
  pheno_tune <- left_join(pheno_tune,rarevariant_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:17],"RV_PRS")
  
  PRSs_Tune <- pheno_tune[,c(17,18)]
  
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tune)
  
  PRSs_Tune$Residuals <- model.null$residuals
  
  model.best <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
  colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt")
  
  rarevariant_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_Burden_Validation_All.txt")
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:16],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:17],"RV_PRS")
  
  PRSs_Validation <- pheno_vad[,c(17,18)]
  
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad)
  
  PRSs_Validation$Residuals <- model.null$residuals
  
  predicted_prs <- predict(model.best,PRSs_Validation)
  
  r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square
  
  effects <- summary(model.best)$coefficients[-1,1]
  
  save(effects,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Coefficients_All_Burden.RData")
  
  data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "CV_plus_RV_Burden",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Burden_All_Result.RData")  
}
