rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

trait <- "BMI"

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  arrayid <- as.numeric(commandArgs(TRUE)[1])
  
  if(arrayid == 1){
    ## STAARO
    
    ## Pull in Phenotypes/Covariates 
    pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
    
    common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"))
    
    rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_Best_All_STAARO_Tune_All.txt"))
    
    pheno_tuning <- left_join(pheno_tuning,common_prs,by = "IID")
    colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:26],"CV_PRS")
    pheno_tuning <- left_join(pheno_tuning,rarevariant_prs,by = "IID")
    colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:27],"RV_PRS")
    
    PRSs_Tune <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(27,28)]
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning)
    
    PRSs_Tune$Residuals <- model.null$residuals
    
    model.full <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)
    
    model.CV <- lm(Residuals~CV_PRS,data = PRSs_Tune)
    
    model.RV <- lm(Residuals~RV_PRS,data = PRSs_Tune)
    
    model.best <- list(model.full,model.CV,model.RV)[[which.max(c(summary(model.full)$r.square,summary(model.CV)$r.square,summary(model.RV)$r.square))]]
    
    pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
    
    common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
    
    rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_Best_All_STAARO_Validation_All.txt"))
    
    pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
    colnames(pheno_vad) <- c(colnames(pheno_vad)[1:26],"CV_PRS")
    pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
    colnames(pheno_vad) <- c(colnames(pheno_vad)[1:27],"RV_PRS")
    
    pheno_vad <- pheno_vad[!is.na(pheno_vad[,trait]),]
    PRSs_Validation <- pheno_vad[,c(1,27,28)]
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    PRSs_Validation_EUR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    PRSs_Validation_NonEur <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    PRSs_Validation_UNK <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    PRSs_Validation_SAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    PRSs_Validation_MIX <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    PRSs_Validation_AFR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    PRSs_Validation_EAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    #EUR
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    
    predicted_prs <- predict(model.best,PRSs_Validation_EUR)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    effects <- summary(model.best)$coefficients[-1,1]
    
    save(effects,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Coefficients_All_STAARO.RData"))
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_EUR",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_EUR.RData"))
    
    #NonEur
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    
    predicted_prs <- predict(model.best,PRSs_Validation_NonEur)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_NonEur",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_NonEur.RData"))
    
    #AFR
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    
    predicted_prs <- predict(model.best,PRSs_Validation_AFR)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_AFR",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_AFR.RData"))
    
    
    #EAS
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    
    predicted_prs <- predict(model.best,PRSs_Validation_EAS)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_EAS",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_EAS.RData"))
    
    #SAS
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    
    predicted_prs <- predict(model.best,PRSs_Validation_SAS)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_SAS",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_SAS.RData"))
    
    #MIX
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    
    predicted_prs <- predict(model.best,PRSs_Validation_MIX)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_MIX",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_MIX.RData"))
    
    
    #UNK
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    
    predicted_prs <- predict(model.best,PRSs_Validation_UNK)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_UNK",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_UNK.RData"))
    
  }else{
    ## Burden
    
    ## Pull in Phenotypes/Covariates 
    pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
    
    common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"))
    
    rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_Best_All_Burden_Tune_All.txt"))
    
    pheno_tuning <- left_join(pheno_tuning,common_prs,by = "IID")
    colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:26],"CV_PRS")
    pheno_tuning <- left_join(pheno_tuning,rarevariant_prs,by = "IID")
    colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:27],"RV_PRS")
    
    PRSs_Tune <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(27,28)]
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning)
    
    PRSs_Tune$Residuals <- model.null$residuals
    
    model.full <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)
    
    model.CV <- lm(Residuals~CV_PRS,data = PRSs_Tune)
    
    model.RV <- lm(Residuals~RV_PRS,data = PRSs_Tune)
    
    model.best <- list(model.full,model.CV,model.RV)[[which.max(c(summary(model.full)$r.square,summary(model.CV)$r.square,summary(model.RV)$r.square))]]
    
    pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
    
    common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
    
    rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_Best_All_Burden_Validation_All.txt"))
    
    pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
    colnames(pheno_vad) <- c(colnames(pheno_vad)[1:26],"CV_PRS")
    pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
    colnames(pheno_vad) <- c(colnames(pheno_vad)[1:27],"RV_PRS")
    
    pheno_vad <- pheno_vad[!is.na(pheno_vad[,trait]),]
    PRSs_Validation <- pheno_vad[,c(1,27,28)]
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    PRSs_Validation_EUR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    PRSs_Validation_NonEur <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    PRSs_Validation_UNK <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    PRSs_Validation_SAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    PRSs_Validation_MIX <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    PRSs_Validation_AFR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    PRSs_Validation_EAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    #EUR
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    
    predicted_prs <- predict(model.best,PRSs_Validation_EUR)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    effects <- summary(model.best)$coefficients[-1,1]
    
    save(effects,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Coefficients_All_Burden.RData"))
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_EUR",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_EUR.RData"))
    
    #NonEur
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    
    predicted_prs <- predict(model.best,PRSs_Validation_NonEur)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_NonEur",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_NonEur.RData"))
    
    #AFR
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    
    predicted_prs <- predict(model.best,PRSs_Validation_AFR)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_AFR",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_AFR.RData"))
    
    
    #EAS
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    
    predicted_prs <- predict(model.best,PRSs_Validation_EAS)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_EAS",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_EAS.RData"))
    
    #SAS
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    
    predicted_prs <- predict(model.best,PRSs_Validation_SAS)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_SAS",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_SAS.RData"))
    
    #MIX
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    
    predicted_prs <- predict(model.best,PRSs_Validation_MIX)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_MIX",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_MIX.RData"))
    
    
    #UNK
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    
    predicted_prs <- predict(model.best,PRSs_Validation_UNK)
    
    r2 <- summary(lm(model.null$residuals~predicted_prs))$r.square
    
    data <- data.frame(y = model.null$residuals, x = predicted_prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    SL.result <- data.frame(method = "CV_plus_RV_Burden_UNK",
                            r2 = r2,
                            r2_low = ci_result$percent[4],
                            r2_high = ci_result$percent[5]
    )
    
    save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_UNK.RData")) 
  }
}