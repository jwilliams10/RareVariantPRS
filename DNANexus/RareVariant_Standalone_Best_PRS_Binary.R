rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/RareVariant_Standalone_Best_PRS_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/RareVariant_Standalone_Best_PRS_Binary.sh -icmd="bash RareVariant_Standalone_Best_PRS_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/ --priority low --instance-type mem1_ssd1_v2_x4

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){

  Noncoding_Burden <- read.csv(paste0(trait,"_Noncoding_Burden_PRS.csv"))
  Noncoding_STAARO <- read.csv(paste0(trait,"_Noncoding_STAARO_PRS.csv"))
  
  Coding_Burden <- read.csv(paste0(trait,"_Coding_Burden_PRS.csv"))
  Coding_STAARO <- read.csv(paste0(trait,"_Coding_STAARO_PRS.csv"))
  
  ncRNA_Burden <- read.csv(paste0(trait,"_ncRNA_Burden_PRS.csv"))
  ncRNA_STAARO <- read.csv(paste0(trait,"_ncRNA_STAARO_PRS.csv"))
  
  Tune_Null_Model <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
  Validation_Null_Model <- get(load(paste0(trait,"_Validation_Null_Model.RData")))
  
  file.remove(paste0(trait,"_Noncoding_Burden_PRS.csv"))
  file.remove(paste0(trait,"_Noncoding_STAARO_PRS.csv"))
  file.remove(paste0(trait,"_Coding_Burden_PRS.csv"))
  file.remove(paste0(trait,"_Coding_STAARO_PRS.csv"))
  file.remove(paste0(trait,"_ncRNA_Burden_PRS.csv"))
  file.remove(paste0(trait,"_ncRNA_STAARO_PRS.csv"))
  file.remove(paste0(trait,"_Tune_Null_Model.RData"))
  file.remove(paste0(trait,"_Validation_Null_Model.RData"))
  
  
  for(i in 2:16){
    Noncoding_Burden[,i] <- Noncoding_Burden[,i] + ncRNA_Burden[,i]
    Noncoding_STAARO[,i] <- Noncoding_STAARO[,i] + ncRNA_STAARO[,i]
  }
  rm(ncRNA_Burden);rm(ncRNA_STAARO)
  
  
  
  Noncoding_Burden_Tune <- Noncoding_Burden[Noncoding_Burden$IID %in% Tune_Null_Model$id_include,]
  Noncoding_Burden_Validation <- Noncoding_Burden[Noncoding_Burden$IID %in% Validation_Null_Model$id_include,]
  
  Noncoding_STAARO_Tune <- Noncoding_STAARO[Noncoding_STAARO$IID %in% Tune_Null_Model$id_include,]
  Noncoding_STAARO_Validation <- Noncoding_STAARO[Noncoding_STAARO$IID %in% Validation_Null_Model$id_include,]
  
  
  Coding_Burden_Tune <- Coding_Burden[Coding_Burden$IID %in% Tune_Null_Model$id_include,]
  Coding_Burden_Validation <- Coding_Burden[Coding_Burden$IID %in% Validation_Null_Model$id_include,]
  
  Coding_STAARO_Tune <- Coding_STAARO[Coding_STAARO$IID %in% Tune_Null_Model$id_include,]
  Coding_STAARO_Validation <- Coding_STAARO[Coding_STAARO$IID %in% Validation_Null_Model$id_include,]
  
  
  
  Burden_Combined_Tune <- cbind(Noncoding_Burden_Tune,Coding_Burden_Tune[,-1])
  colnames(Burden_Combined_Tune) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  STAARO_Combined_Tune <- cbind(Noncoding_STAARO_Tune,Coding_STAARO_Tune[,-1])
  colnames(STAARO_Combined_Tune) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  
  Burden_Combined_Validation <- cbind(Noncoding_Burden_Validation,Coding_Burden_Validation[,-1])
  colnames(Burden_Combined_Validation) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  STAARO_Combined_Validation <- cbind(Noncoding_STAARO_Validation,Coding_STAARO_Validation[,-1])
  colnames(STAARO_Combined_Validation) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  
  
  
  ## Pull in Phenotypes/Covariates 
  pheno_tuning <- read.delim("All_Tune.txt")
  
  pheno_tuning_STAARO <- left_join(pheno_tuning,STAARO_Combined_Tune,by = "IID")
  pheno_tuning_STAARO <- pheno_tuning_STAARO[!is.na(pheno_tuning_STAARO[,trait]),]
  STAARO_Combined_Tune <- pheno_tuning_STAARO[,-c(1:26)]
  
  pheno_tuning_Burden <- left_join(pheno_tuning,Burden_Combined_Tune,by = "IID")
  pheno_tuning_Burden <- pheno_tuning_Burden[!is.na(pheno_tuning_Burden[,trait]),]
  Burden_Combined_Tune <- pheno_tuning_Burden[,-c(1:26)]
  
  pheno_vad <- read.delim("All_Validation.txt")
  
  pheno_vad_STAARO <- left_join(pheno_vad,STAARO_Combined_Validation,by = "IID")
  pheno_vad_STAARO <- pheno_vad_STAARO[!is.na(pheno_vad_STAARO[,trait]),]
  STAARO_Combined_Validation <- pheno_vad_STAARO[,-c(1:26)]
  
  pheno_vad_Burden <- left_join(pheno_vad,Burden_Combined_Validation,by = "IID")
  pheno_vad_Burden <- pheno_vad_Burden[!is.na(pheno_vad_Burden[,trait]),]
  Burden_Combined_Validation <- pheno_vad_Burden[,-c(1:26)]
  
  
  
  
  
  load("all_phenotypes.RData")
  
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
                          variable =paste0("Coding_Threshold",k),
                          confounders = confounders,
                          data = pheno_tuning_STAARO,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_tun_GeneCentric_Coding_STAARO[k] <- roc_obj$auc
  }
  
  ### EUR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_EUR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EUR",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_EUR_result.RData"))
  
  ### NonEur
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_NonEur,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_NonEur",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_NonEur_result.RData"))
  
  ### EAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_EAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EAS",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_EAS_result.RData"))
  
  ### AFR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_AFR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_AFR",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_AFR_result.RData"))
  
  ### SAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_SAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_SAS",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_SAS_result.RData"))
  
  ### MIX
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_MIX,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_MIX",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_MIX_result.RData"))
  
  ### UNK
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_UNK,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_UNK",
                                                     AUC = AUC_GeneCentric_Coding_STAARO,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_UNK_result.RData"))
  
  
  
  
  
  
  ##### Gene Centric Noncoding: STAARO
  AUC_tun_GeneCentric_Noncoding_STAARO  <- rep(0,15)
  for(k in 1:15){
    roc_obj <- roc.binary(status = trait,
                          variable =paste0("Noncoding_Threshold",k),
                          confounders = confounders,
                          data = pheno_tuning_STAARO,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_tun_GeneCentric_Noncoding_STAARO[k] <- roc_obj$auc
  }
  
  ### EUR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_EUR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EUR",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_EUR_result.RData"))
  
  ### NonEur
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_NonEur,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_NonEur",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_NonEur_result.RData"))
  
  ### EAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_EAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EAS",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_EAS_result.RData"))
  
  ### AFR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_AFR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_AFR",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_AFR_result.RData"))
  
  ### SAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_SAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_SAS",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_SAS_result.RData"))
  
  ### MIX
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_MIX,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_MIX",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_MIX_result.RData"))
  
  ### UNK
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                        confounders = confounders,
                        data = pheno_vad_UNK,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_STAARO <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_STAARO)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_UNK",
                                                        AUC = AUC_GeneCentric_Noncoding_STAARO,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_UNK_result.RData"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###############################
  ######################Burden
  
  load("all_phenotypes.RData")
  
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
                          variable =paste0("Coding_Threshold",k),
                          confounders = confounders,
                          data = pheno_tuning_Burden,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_tun_GeneCentric_Coding_Burden[k] <- roc_obj$auc
  }
  
  ### EUR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_EUR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EUR",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_EUR_result.RData"))
  
  ### NonEur
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_NonEur,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_NonEur",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_NonEur_result.RData"))
  
  ### EAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_EAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EAS",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_EAS_result.RData"))
  
  ### AFR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_AFR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_AFR",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_AFR_result.RData"))
  
  ### SAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_SAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_SAS",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_SAS_result.RData"))
  
  ### MIX
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_MIX,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_MIX",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_MIX_result.RData"))
  
  ### UNK
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_UNK,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Coding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Coding_Threshold",which.max(AUC_tun_GeneCentric_Coding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_UNK",
                                                     AUC = AUC_GeneCentric_Coding_Burden,
                                                     AUC_low = ci_result$percent[4],
                                                     AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_UNK_result.RData"))
  
  
  
  
  
  
  ##### Gene Centric Noncoding: Burden
  AUC_tun_GeneCentric_Noncoding_Burden  <- rep(0,15)
  for(k in 1:15){
    roc_obj <- roc.binary(status = trait,
                          variable =paste0("Noncoding_Threshold",k),
                          confounders = confounders,
                          data = pheno_tuning_Burden,
                          precision=seq(0.05,0.95, by=0.05))
    AUC_tun_GeneCentric_Noncoding_Burden[k] <- roc_obj$auc
  }
  
  ### EUR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_EUR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EUR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EUR",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_EUR_result.RData"))
  
  ### NonEur
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_NonEur,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_NonEur, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_NonEur",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_NonEur_result.RData"))
  
  ### EAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_EAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_EAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EAS",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_EAS_result.RData"))
  
  ### AFR
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_AFR,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_AFR, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_AFR",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_AFR_result.RData"))
  
  ### SAS
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_SAS,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_SAS, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_SAS",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_SAS_result.RData"))
  
  ### MIX
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_MIX,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_MIX, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_MIX",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_MIX_result.RData"))
  
  ### UNK
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                        confounders = confounders,
                        data = pheno_vad_UNK,
                        precision=seq(0.05,0.95, by=0.05))
  AUC_GeneCentric_Noncoding_Burden <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("Noncoding_Threshold",which.max(AUC_tun_GeneCentric_Noncoding_Burden)),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = pheno_vad_UNK, statistic = calc_auc, R = 1000)
  
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  AUC.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_UNK",
                                                        AUC = AUC_GeneCentric_Noncoding_Burden,
                                                        AUC_low = ci_result$percent[4],
                                                        AUC_high = ci_result$percent[5]
  )
  
  save(AUC.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_UNK_result.RData"))
  
  
}

system("rm All_Tune.txt")
system("rm All_Validation.txt")
system("rm all_phenotypes.RData")