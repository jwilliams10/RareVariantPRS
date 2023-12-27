rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

trait <- "BMI"

for(trait in c("BMI","HDL","LDL","Height","TC","logTG")){
  STAARO_GeneCentric_Coding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_STAARO_GeneCentric_Coding_Tune_PRS.csv"))
  STAARO_GeneCentric_Coding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_STAARO_GeneCentric_Coding_Validation_PRS.csv"))
  
  STAARO_GeneCentric_Noncoding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_STAARO_GeneCentric_Noncoding_Tune_PRS.csv"))
  STAARO_GeneCentric_Noncoding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_STAARO_GeneCentric_Noncoding_Validation_PRS.csv"))
  
  STAARO_SlidingWindow_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_STAARO_SlidingWindow_Tune_PRS.csv"))
  STAARO_SlidingWindow_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_STAARO_SlidingWindow_Validation_PRS.csv"))
  
  colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Coding_Tune_PRS)]))
  colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Tune_PRS)]))
  colnames(STAARO_SlidingWindow_Tune_PRS) <- c(colnames(STAARO_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Tune_PRS)[2:ncol(STAARO_SlidingWindow_Tune_PRS)]))
  STAARO_Combined_Tune <- cbind(STAARO_GeneCentric_Coding_Tune_PRS,STAARO_GeneCentric_Noncoding_Tune_PRS[,-1],STAARO_SlidingWindow_Tune_PRS[,-1])
  
  colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Coding_Validation_PRS)]))
  colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Validation_PRS)]))
  colnames(STAARO_SlidingWindow_Validation_PRS) <- c(colnames(STAARO_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Validation_PRS)[2:ncol(STAARO_SlidingWindow_Validation_PRS)]))
  STAARO_Combined_Validation <- cbind(STAARO_GeneCentric_Coding_Validation_PRS,STAARO_GeneCentric_Noncoding_Validation_PRS[,-1],STAARO_SlidingWindow_Validation_PRS[,-1])
  
  
  Burden_GeneCentric_Coding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_Burden_GeneCentric_Coding_Tune_PRS.csv"))
  Burden_GeneCentric_Coding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_Burden_GeneCentric_Coding_Validation_PRS.csv"))
  
  Burden_GeneCentric_Noncoding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_Burden_GeneCentric_Noncoding_Tune_PRS.csv"))
  Burden_GeneCentric_Noncoding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_Burden_GeneCentric_Noncoding_Validation_PRS.csv"))
  
  Burden_SlidingWindow_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_Burden_SlidingWindow_Tune_PRS.csv"))
  Burden_SlidingWindow_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_Burden_SlidingWindow_Validation_PRS.csv"))
  
  colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Tune_PRS)[2:ncol(Burden_GeneCentric_Coding_Tune_PRS)]))
  colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Tune_PRS)]))
  colnames(Burden_SlidingWindow_Tune_PRS) <- c(colnames(Burden_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Tune_PRS)[2:ncol(Burden_SlidingWindow_Tune_PRS)]))
  Burden_Combined_Tune <- cbind(Burden_GeneCentric_Coding_Tune_PRS,Burden_GeneCentric_Noncoding_Tune_PRS[,-1],Burden_SlidingWindow_Tune_PRS[,-1])
  
  colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Validation_PRS)[2:ncol(Burden_GeneCentric_Coding_Validation_PRS)]))
  colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Validation_PRS)]))
  colnames(Burden_SlidingWindow_Validation_PRS) <- c(colnames(Burden_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Validation_PRS)[2:ncol(Burden_SlidingWindow_Validation_PRS)]))
  Burden_Combined_Validation <- cbind(Burden_GeneCentric_Coding_Validation_PRS,Burden_GeneCentric_Noncoding_Validation_PRS[,-1],Burden_SlidingWindow_Validation_PRS[,-1])
  
  ## Pull in Phenotypes/Covariates 
  pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  
  pheno_tuning_STAARO <- left_join(pheno_tuning,STAARO_Combined_Tune,by = "IID")
  pheno_tuning_STAARO <- pheno_tuning_STAARO[!is.na(pheno_tuning_STAARO[,trait]),]
  STAARO_Combined_Tune <- pheno_tuning_STAARO[,-c(1:26)]
  
  pheno_tuning_Burden <- left_join(pheno_tuning,Burden_Combined_Tune,by = "IID")
  pheno_tuning_Burden <- pheno_tuning_Burden[!is.na(pheno_tuning_Burden[,trait]),]
  Burden_Combined_Tune <- pheno_tuning_Burden[,-c(1:26)]
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  
  pheno_vad_STAARO <- left_join(pheno_vad,STAARO_Combined_Validation,by = "IID")
  pheno_vad_STAARO <- pheno_vad_STAARO[!is.na(pheno_vad_STAARO[,trait]),]
  STAARO_Combined_Validation <- pheno_vad_STAARO[,-c(1:26)]
  
  pheno_vad_Burden <- left_join(pheno_vad,Burden_Combined_Validation,by = "IID")
  pheno_vad_Burden <- pheno_vad_Burden[!is.na(pheno_vad_Burden[,trait]),]
  Burden_Combined_Validation <- pheno_vad_Burden[,-c(1:26)]
  
  ############### R2's
  arrayid <- as.numeric(commandArgs(TRUE)[1])
  
  if(arrayid == 1){
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_vad_EUR <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_vad_NonEur <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    pheno_vad_UNK <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    pheno_vad_SAS <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_vad_MIX <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    pheno_vad_AFR <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_vad_EAS <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    ##### Gene Centric Coding: STAARO
    r2_tun_GeneCentric_Coding_STAARO  <- rep(0,15)
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_STAARO)
    for(k in 1:15){
      prs <- pheno_tuning_STAARO[,paste0("GeneCentric_Coding_PRS_Threshold_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning_STAARO)
      r2_tun_GeneCentric_Coding_STAARO[k] <- summary(model.prs)$r.square
    }
    
    ### EUR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    prs <- pheno_vad_EUR[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EUR",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EUR_result.RData"))
    
    ### NonEur
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    prs <- pheno_vad_NonEur[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_NonEur",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_NonEur_result.RData"))
    
    ### EAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    prs <- pheno_vad_EAS[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EAS",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EAS_result.RData"))
    
    ### AFR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    prs <- pheno_vad_AFR[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_AFR",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_AFR_result.RData"))
    
    ### SAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    prs <- pheno_vad_SAS[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_SAS",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_SAS_result.RData"))
    
    ### MIX
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    prs <- pheno_vad_MIX[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_MIX",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_MIX_result.RData"))
    
    ### UNK
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    prs <- pheno_vad_UNK[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_UNK",
                                                      r2 = r2_GeneCentric_Coding_STAARO,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_UNK_result.RData"))
    
    
    
    
    
    
    ##### Gene Centric Noncoding: STAARO
    r2_tun_GeneCentric_Noncoding_STAARO  <- rep(0,15)
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_STAARO)
    for(k in 1:15){
      prs <- pheno_tuning_STAARO[,paste0("GeneCentric_Noncoding_PRS_Threshold_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning_STAARO)
      r2_tun_GeneCentric_Noncoding_STAARO[k] <- summary(model.prs)$r.square
    }
    
    ### EUR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    prs <- pheno_vad_EUR[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EUR",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EUR_result.RData"))
    
    ### NonEur
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    prs <- pheno_vad_NonEur[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_NonEur",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_NonEur_result.RData"))
    
    ### EAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    prs <- pheno_vad_EAS[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EAS",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EAS_result.RData"))
    
    ### AFR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    prs <- pheno_vad_AFR[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_AFR",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_AFR_result.RData"))
    
    ### SAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    prs <- pheno_vad_SAS[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_SAS",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_SAS_result.RData"))
    
    ### MIX
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    prs <- pheno_vad_MIX[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_MIX",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_MIX_result.RData"))
    
    ### UNK
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    prs <- pheno_vad_UNK[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_UNK",
                                                         r2 = r2_GeneCentric_Noncoding_STAARO,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_UNK_result.RData"))
    
    ##### Sliding Window: STAARO
    r2_tun_SlidingWindow_STAARO  <- rep(0,15)
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_STAARO)
    for(k in 1:15){
      prs <- pheno_tuning_STAARO[,paste0("SlidingWindow_PRS_Threshold_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning_STAARO)
      r2_tun_SlidingWindow_STAARO[k] <- summary(model.prs)$r.square
    }
    
    ### EUR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    prs <- pheno_vad_EUR[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_EUR",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EUR_result.RData"))
    
    ### NonEur
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    prs <- pheno_vad_NonEur[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_NonEur",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_NonEur_result.RData"))
    
    ### EAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    prs <- pheno_vad_EAS[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_EAS",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EAS_result.RData"))
    
    ### AFR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    prs <- pheno_vad_AFR[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_AFR",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_AFR_result.RData"))
    
    ### SAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    prs <- pheno_vad_SAS[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_SAS",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_SAS_result.RData"))
    
    ### MIX
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    prs <- pheno_vad_MIX[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_MIX",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_MIX_result.RData"))
    
    ### UNK
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    prs <- pheno_vad_UNK[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_UNK",
                                                 r2 = r2_SlidingWindow_STAARO,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_UNK_result.RData"))
    
  }else{
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_vad_EUR <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_vad_NonEur <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    pheno_vad_UNK <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    pheno_vad_SAS <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_vad_MIX <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    pheno_vad_AFR <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_vad_EAS <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    ##### Gene Centric Coding: Burden
    r2_tun_GeneCentric_Coding_Burden  <- rep(0,15)
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_Burden)
    for(k in 1:15){
      prs <- pheno_tuning_Burden[,paste0("GeneCentric_Coding_PRS_Threshold_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning_Burden)
      r2_tun_GeneCentric_Coding_Burden[k] <- summary(model.prs)$r.square
    }
    
    ### EUR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    prs <- pheno_vad_EUR[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EUR",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EUR_result.RData"))
    
    ### NonEur
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    prs <- pheno_vad_NonEur[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_NonEur",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_NonEur_result.RData"))
    
    ### EAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    prs <- pheno_vad_EAS[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EAS",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EAS_result.RData"))
    
    ### AFR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    prs <- pheno_vad_AFR[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_AFR",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_AFR_result.RData"))
    
    ### SAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    prs <- pheno_vad_SAS[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_SAS",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_SAS_result.RData"))
    
    ### MIX
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    prs <- pheno_vad_MIX[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_MIX",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_MIX_result.RData"))
    
    ### UNK
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    prs <- pheno_vad_UNK[,paste0("GeneCentric_Coding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_UNK",
                                                      r2 = r2_GeneCentric_Coding_Burden,
                                                      r2_low = ci_result$percent[4],
                                                      r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_UNK_result.RData"))
    
    ##### Gene Centric Noncoding: Burden
    r2_tun_GeneCentric_Noncoding_Burden  <- rep(0,15)
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_Burden)
    for(k in 1:15){
      prs <- pheno_tuning_Burden[,paste0("GeneCentric_Noncoding_PRS_Threshold_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning_Burden)
      r2_tun_GeneCentric_Noncoding_Burden[k] <- summary(model.prs)$r.square
    }
    
    ### EUR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    prs <- pheno_vad_EUR[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EUR",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EUR_result.RData"))
    
    ### NonEur
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    prs <- pheno_vad_NonEur[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_NonEur",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_NonEur_result.RData"))
    
    ### EAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    prs <- pheno_vad_EAS[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EAS",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EAS_result.RData"))
    
    ### AFR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    prs <- pheno_vad_AFR[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_AFR",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_AFR_result.RData"))
    
    ### SAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    prs <- pheno_vad_SAS[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_SAS",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_SAS_result.RData"))
    
    ### MIX
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    prs <- pheno_vad_MIX[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_MIX",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_MIX_result.RData"))
    
    ### UNK
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    prs <- pheno_vad_UNK[,paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_UNK",
                                                         r2 = r2_GeneCentric_Noncoding_Burden,
                                                         r2_low = ci_result$percent[4],
                                                         r2_high = ci_result$percent[5]
    )
    
    save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_UNK_result.RData"))
    
    ##### Sliding Window: Burden
    r2_tun_SlidingWindow_Burden  <- rep(0,15)
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_Burden)
    for(k in 1:15){
      prs <- pheno_tuning_Burden[,paste0("SlidingWindow_PRS_Threshold_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning_Burden)
      r2_tun_SlidingWindow_Burden[k] <- summary(model.prs)$r.square
    }
    
    ### EUR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
    prs <- pheno_vad_EUR[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_EUR",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EUR_result.RData"))
    
    ### NonEur
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
    prs <- pheno_vad_NonEur[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_NonEur",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_NonEur_result.RData"))
    
    ### EAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
    prs <- pheno_vad_EAS[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_EAS",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EAS_result.RData"))
    
    ### AFR
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
    prs <- pheno_vad_AFR[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_AFR",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_AFR_result.RData"))
    
    ### SAS
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
    prs <- pheno_vad_SAS[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_SAS",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_SAS_result.RData"))
    
    ### MIX
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
    prs <- pheno_vad_MIX[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_MIX",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_MIX_result.RData"))
    
    ### UNK
    model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
    prs <- pheno_vad_UNK[,paste0("SlidingWindow_PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
    model.vad.prs <- lm(model.vad.null$residual~prs)
    r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
    
    data <- data.frame(y = model.vad.null$residual, x = prs)
    R2Boot <- function(data,indices){
      boot_data <- data[indices, ]
      model <- lm(y ~ x, data = boot_data)
      result <- summary(model)$r.square
      return(c(result))
    }
    library(boot)
    boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
    
    ci_result <- boot.ci(boot_r2, type = "perc")
    r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_UNK",
                                                 r2 = r2_SlidingWindow_Burden,
                                                 r2_low = ci_result$percent[4],
                                                 r2_high = ci_result$percent[5]
    )
    
    save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_UNK_result.RData"))
  }
}



