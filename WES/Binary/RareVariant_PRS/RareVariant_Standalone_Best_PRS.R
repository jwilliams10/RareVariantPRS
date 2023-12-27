rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  STAARO_GeneCentric_Coding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_STAARO_GeneCentric_Coding_Tune_PRS.csv"))
  STAARO_GeneCentric_Coding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_STAARO_GeneCentric_Coding_Validation_PRS.csv"))
  
  STAARO_GeneCentric_Noncoding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_STAARO_GeneCentric_Noncoding_Tune_PRS.csv"))
  STAARO_GeneCentric_Noncoding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_STAARO_GeneCentric_Noncoding_Validation_PRS.csv"))
  
  STAARO_SlidingWindow_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_STAARO_SlidingWindow_Tune_PRS.csv"))
  STAARO_SlidingWindow_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_STAARO_SlidingWindow_Validation_PRS.csv"))
  
  colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Coding_Tune_PRS)]))
  colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Tune_PRS)]))
  colnames(STAARO_SlidingWindow_Tune_PRS) <- c(colnames(STAARO_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Tune_PRS)[2:ncol(STAARO_SlidingWindow_Tune_PRS)]))
  STAARO_Combined_Tune <- cbind(STAARO_GeneCentric_Coding_Tune_PRS,STAARO_GeneCentric_Noncoding_Tune_PRS[,-1],STAARO_SlidingWindow_Tune_PRS[,-1])
  
  colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Coding_Validation_PRS)]))
  colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Validation_PRS)]))
  colnames(STAARO_SlidingWindow_Validation_PRS) <- c(colnames(STAARO_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Validation_PRS)[2:ncol(STAARO_SlidingWindow_Validation_PRS)]))
  STAARO_Combined_Validation <- cbind(STAARO_GeneCentric_Coding_Validation_PRS,STAARO_GeneCentric_Noncoding_Validation_PRS[,-1],STAARO_SlidingWindow_Validation_PRS[,-1])
  
  
  Burden_GeneCentric_Coding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_Burden_GeneCentric_Coding_Tune_PRS.csv"))
  Burden_GeneCentric_Coding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_Burden_GeneCentric_Coding_Validation_PRS.csv"))
  
  Burden_GeneCentric_Noncoding_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_Burden_GeneCentric_Noncoding_Tune_PRS.csv"))
  Burden_GeneCentric_Noncoding_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_Burden_GeneCentric_Noncoding_Validation_PRS.csv"))
  
  Burden_SlidingWindow_Tune_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_Burden_SlidingWindow_Tune_PRS.csv"))
  Burden_SlidingWindow_Validation_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_Burden_SlidingWindow_Validation_PRS.csv"))
  
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
  
  ############### AUC's
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
    
    if(trait %in% c("Breast","Prostate")){
      confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    }else{
      confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    }
    
    ##### Gene Centric Coding: STAARO
    AUC_tun_GeneCentric_Coding_STAARO  <- rep(0,15)
    for(k in 1:15){
      roc_obj <- roc.binary(status = trait,
                            variable =paste0("GeneCentric_Coding_PRS_Threshold_",k),
                            confounders = confounders,
                            data = pheno_tuning_STAARO,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_tun_GeneCentric_Coding_STAARO[k] <- roc_obj$auc
    }
    
    ### EUR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_EUR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EUR",
                                                      AUC = AUC_GeneCentric_Coding_STAARO,
                                                      AUC_low = ci_result$percent[4],
                                                      AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EUR_result.RData"))
    
    ### NonEur
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_NonEur,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_NonEur",
                                                       AUC = AUC_GeneCentric_Coding_STAARO,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_NonEur_result.RData"))
    
    ### EAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_EAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EAS",
                                                       AUC = AUC_GeneCentric_Coding_STAARO,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EAS_result.RData"))
    
    ### AFR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_AFR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_AFR",
                                                       AUC = AUC_GeneCentric_Coding_STAARO,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_AFR_result.RData"))
    
    ### SAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_SAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_SAS",
                                                       AUC = AUC_GeneCentric_Coding_STAARO,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_SAS_result.RData"))
    
    ### MIX
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_MIX,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_MIX",
                                                       AUC = AUC_GeneCentric_Coding_STAARO,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_MIX_result.RData"))
    
    ### UNK
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_UNK,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_UNK",
                                                       AUC = AUC_GeneCentric_Coding_STAARO,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_UNK_result.RData"))
    
    
    
    
    
    
    ##### Gene Centric Noncoding: STAARO
    AUC_tun_GeneCentric_Noncoding_STAARO  <- rep(0,15)
    for(k in 1:15){
      roc_obj <- roc.binary(status = trait,
                            variable =paste0("GeneCentric_Noncoding_PRS_Threshold_",k),
                            confounders = confounders,
                            data = pheno_tuning_STAARO,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_tun_GeneCentric_Noncoding_STAARO[k] <- roc_obj$auc
    }
    
    ### EUR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_EUR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EUR",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EUR_result.RData"))
    
    ### NonEur
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_NonEur,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_NonEur",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_NonEur_result.RData"))
    
    ### EAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_EAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EAS",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EAS_result.RData"))
    
    ### AFR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_AFR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_AFR",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_AFR_result.RData"))
    
    ### SAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_SAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_SAS",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_SAS_result.RData"))
    
    ### MIX
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_MIX,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_MIX",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_MIX_result.RData"))
    
    ### UNK
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_UNK,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_UNK",
                                                          AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_UNK_result.RData"))
    
    
    ##### Sliding Window: STAARO
    AUC_tun_SlidingWindow_STAARO  <- rep(0,15)
    for(k in 1:15){
      roc_obj <- roc.binary(status = trait,
                            variable =paste0("SlidingWindow_PRS_Threshold_",k),
                            confounders = confounders,
                            data = pheno_tuning_STAARO,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_tun_SlidingWindow_STAARO[k] <- roc_obj$auc
    }
    
    ### EUR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_EUR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_EUR",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EUR_result.RData"))
    
    ### NonEur
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_NonEur,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_NonEur",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_NonEur_result.RData"))
    
    ### EAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_EAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_EAS",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EAS_result.RData"))
    
    ### AFR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_AFR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_AFR",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_AFR_result.RData"))
    
    ### SAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_SAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_SAS",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_SAS_result.RData"))
    
    ### MIX
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_MIX,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_MIX",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_MIX_result.RData"))
    
    ### UNK
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                          confounders = confounders,
                          data = pheno_vad_UNK,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_STAARO <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_STAARO)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_UNK",
                                                  AUC = AUC_SlidingWindow_STAARO,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_UNK_result.RData"))
    
  }else{
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_vad_EUR <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_vad_NonEur <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
    pheno_vad_UNK <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
    pheno_vad_SAS <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_vad_MIX <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
    pheno_vad_AFR <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_vad_EAS <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    if(trait %in% c("Breast","Prostate")){
      confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    }else{
      confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    }
    
    ##### Gene Centric Coding: Burden
    AUC_tun_GeneCentric_Coding_Burden  <- rep(0,15)
    for(k in 1:15){
      roc_obj <- roc.binary(status = trait,
                            variable =paste0("GeneCentric_Coding_PRS_Threshold_",k),
                            confounders = confounders,
                            data = pheno_tuning_Burden,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_tun_GeneCentric_Coding_Burden[k] <- roc_obj$auc
    }
    
    ### EUR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_EUR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EUR",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EUR_result.RData"))
    
    ### NonEur
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_NonEur,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_NonEur",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_NonEur_result.RData"))
    
    ### EAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_EAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EAS",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EAS_result.RData"))
    
    ### AFR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_AFR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_AFR",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_AFR_result.RData"))
    
    ### SAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_SAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_SAS",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_SAS_result.RData"))
    
    ### MIX
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_MIX,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_MIX",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_MIX_result.RData"))
    
    ### UNK
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_UNK,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Coding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Coding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_UNK",
                                                       AUC = AUC_GeneCentric_Coding_Burden,
                                                       AUC_low = ci_result$percent[4],
                                                       AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_UNK_result.RData"))
    
    
    
    
    
    
    ##### Gene Centric Noncoding: Burden
    AUC_tun_GeneCentric_Noncoding_Burden  <- rep(0,15)
    for(k in 1:15){
      roc_obj <- roc.binary(status = trait,
                            variable =paste0("GeneCentric_Noncoding_PRS_Threshold_",k),
                            confounders = confounders,
                            data = pheno_tuning_Burden,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_tun_GeneCentric_Noncoding_Burden[k] <- roc_obj$auc
    }
    
    ### EUR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_EUR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EUR",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EUR_result.RData"))
    
    ### NonEur
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_NonEur,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_NonEur",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_NonEur_result.RData"))
    
    ### EAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_EAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EAS",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EAS_result.RData"))
    
    ### AFR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_AFR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_AFR",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_AFR_result.RData"))
    
    ### SAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_SAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_SAS",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_SAS_result.RData"))
    
    ### MIX
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_MIX,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_MIX",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_MIX_result.RData"))
    
    ### UNK
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = pheno_vad_UNK,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("GeneCentric_Noncoding_PRS_Threshold_",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_UNK",
                                                          AUC = AUC_GeneCentric_Noncoding_Burden,
                                                          AUC_low = ci_result$percent[4],
                                                          AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_UNK_result.RData"))
    
    
    ##### Sliding Window: Burden
    AUC_tun_SlidingWindow_Burden  <- rep(0,15)
    for(k in 1:15){
      roc_obj <- roc.binary(status = trait,
                            variable =paste0("SlidingWindow_PRS_Threshold_",k),
                            confounders = confounders,
                            data = pheno_tuning_Burden,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_tun_SlidingWindow_Burden[k] <- roc_obj$auc
    }
    
    ### EUR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_EUR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_EUR",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EUR_result.RData"))
    
    ### NonEur
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_NonEur,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_NonEur",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_NonEur_result.RData"))
    
    ### EAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_EAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_EAS",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EAS_result.RData"))
    
    ### AFR
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_AFR,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_AFR",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_AFR_result.RData"))
    
    ### SAS
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_SAS,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_SAS",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_SAS_result.RData"))
    
    ### MIX
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_MIX,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_MIX",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_MIX_result.RData"))
    
    ### UNK
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                          confounders = confounders,
                          data = pheno_vad_UNK,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_SlidingWindow_Burden <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("SlidingWindow_PRS_Threshold_",which.max(AUC_tun_SlidingWindow_Burden)),
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
    
    ci_result <- boot.ci(boot_AUC, type = "perc")
    AUC.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_UNK",
                                                  AUC = AUC_SlidingWindow_Burden,
                                                  AUC_low = ci_result$percent[4],
                                                  AUC_high = ci_result$percent[5]
    )
    
    save(AUC.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_UNK_result.RData"))
  }
}



