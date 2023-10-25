rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- 1

for(i in 1:length(Y_tune)){
  ## STAARO
  
  pheno_tune <- Y_tune[[i]]
  colnames(pheno_tune) <- c("IID","Y")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"))
  
  rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_STAARO_Tune_All",i,".txt"))
  
  pheno_tune <- left_join(pheno_tune,common_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:2],"CV_PRS")
  pheno_tune <- left_join(pheno_tune,rarevariant_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:3],"RV_PRS")
  
  PRSs_Tune <- pheno_tune[,c(3:4)]
  
  model.null <- lm(Y~1,data=pheno_tune)
  
  PRSs_Tune$Residuals <- model.null$residuals
  
  model.best <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)
  
  load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
  pheno_vad <- Y_validation[[i]]
  colnames(pheno_vad) <- c("IID","Y")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))
  
  rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_STAARO_Validation_All",i,".txt"))
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:2],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:3],"RV_PRS")
  
  PRSs_Validation <- pheno_vad[,c(3,4)]
  
  model.null <- lm(Y~1,data=pheno_vad)
  
  PRSs_Validation$Residuals <- model.null$residuals
  
  predicted_prs <- predict(model.best,PRSs_Validation)
  
  r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square
  
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
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result",i,".RData"))
  
  ## Burden
  
  pheno_tune <- Y_tune[[i]]
  colnames(pheno_tune) <- c("IID","Y")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"))
  
  rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_Burden_Tune_All",i,".txt"))
  
  pheno_tune <- left_join(pheno_tune,common_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:2],"CV_PRS")
  pheno_tune <- left_join(pheno_tune,rarevariant_prs,by = "IID")
  colnames(pheno_tune) <- c(colnames(pheno_tune)[1:3],"RV_PRS")
  
  PRSs_Tune <- pheno_tune[,c(3,4)]
  
  model.null <- lm(Y~1,data=pheno_tune)
  
  PRSs_Tune$Residuals <- model.null$residuals
  
  model.best <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)
  
  load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
  pheno_vad <- Y_validation[[i]]
  colnames(pheno_vad) <- c("IID","Y")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))
  
  rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_Burden_Validation_All",i,".txt"))
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:2],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:3],"RV_PRS")
  
  PRSs_Validation <- pheno_vad[,c(3,4)]
  
  model.null <- lm(Y~1,data=pheno_vad)
  
  PRSs_Validation$Residuals <- model.null$residuals
  
  predicted_prs <- predict(model.best,PRSs_Validation)
  
  r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square
  
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
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result",i,".RData"))
  
}