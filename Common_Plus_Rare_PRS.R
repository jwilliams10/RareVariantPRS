rm(list = ls())

arrayid <- as.numeric(commandArgs(TRUE)[1])

if(arrayid == 1){
  ## STAARO
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
  colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt")
  
  rarevariant_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_STAARO_Validation_All.txt")
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:16],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:17],"RV_PRS")
  
  PRSs <- pheno_vad[,c(17,18)]
  
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad)
  
  PRSs$Residuals <- model.null$residuals
  
  r2 <- summary(lm(Residuals~CV_PRS + RV_PRS,data = PRSs))$r.square
  
  effects <- summary(lm(Residuals~CV_PRS + RV_PRS,data = PRSs))$coefficients[-1,1]
  
  save(effects,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Coefficients.RData")
  
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ CV_PRS + RV_PRS, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = PRSs, statistic = R2Boot, R = 100000)
  
  ci_result <- boot.ci(boot_r2, type = "bca")
  SL.result <- data.frame(method = "CV_plus_RV_STAARO",
                          r2 = r2,
                          r2_low = ci_result$bca[4],
                          r2_high = ci_result$bca[5]
  )
  
  save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/STAARO_Result.RData") 
}else{
  ## Burden
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
  colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  
  common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt")
  
  rarevariant_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_Burden_Validation_All.txt")
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:16],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:17],"RV_PRS")
  
  PRSs <- pheno_vad[,c(17,18)]
  
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad)
  
  PRSs$Residuals <- model.null$residuals
  
  r2 <- summary(lm(Residuals~CV_PRS + RV_PRS,data = PRSs))$r.square
  
  effects <- summary(lm(Residuals~CV_PRS + RV_PRS,data = PRSs))$coefficients[-1,1]
  
  save(effects,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Coefficients.RData")
  
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ CV_PRS + RV_PRS, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = PRSs, statistic = R2Boot, R = 100000)
  
  ci_result <- boot.ci(boot_r2, type = "bca")
  SL.result <- data.frame(method = "CV_plus_RV_Burden",
                          r2 = r2,
                          r2_low = ci_result$bca[4],
                          r2_high = ci_result$bca[5]
  )
  
  save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Burden_Result.RData") 
}
