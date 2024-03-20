rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

trait <- "BMI"

for(trait in c("BMI","HDL","LDL","Height","TC","logTG")){

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
  
  
  #########################
  ################### STAARO
  
                                                  
  ##### Gene Centric Coding: STAARO
  r2_tun_GeneCentric_Coding_STAARO  <- rep(0,15)                         
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_STAARO)
  for(k in 1:15){
    prs <- pheno_tuning_STAARO[,paste0("Coding_Threshold",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning_STAARO)
    r2_tun_GeneCentric_Coding_STAARO[k] <- summary(model.prs)$r.square
  }
  
  ### EUR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- pheno_vad_EUR[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_EUR_result.RData"))
  
  ### NonEur
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- pheno_vad_NonEur[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_NonEur_result.RData"))
  
  ### EAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- pheno_vad_EAS[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_EAS_result.RData"))
  
  ### AFR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- pheno_vad_AFR[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_AFR_result.RData"))
  
  ### SAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- pheno_vad_SAS[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_SAS_result.RData"))
  
  ### MIX
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- pheno_vad_MIX[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_MIX_result.RData"))
  
  ### UNK
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- pheno_vad_UNK[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  
  save(r2.result_GeneCentric_Coding_STAARO, file = paste0(trait,"_GeneCentric_Coding_STAARO_UNK_result.RData"))
  
  ##### Gene Centric Noncoding: STAARO
  r2_tun_GeneCentric_Noncoding_STAARO  <- rep(0,15)
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_STAARO)
  for(k in 1:15){
    prs <- pheno_tuning_STAARO[,paste0("Noncoding_Threshold",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning_STAARO)
    r2_tun_GeneCentric_Noncoding_STAARO[k] <- summary(model.prs)$r.square
  }
  
  ### EUR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- pheno_vad_EUR[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_EUR_result.RData"))
  
  ### NonEur
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- pheno_vad_NonEur[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_NonEur_result.RData"))
  
  ### EAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- pheno_vad_EAS[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_EAS_result.RData"))
  
  ### AFR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- pheno_vad_AFR[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_AFR_result.RData"))
  
  ### SAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- pheno_vad_SAS[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_SAS_result.RData"))
  
  ### MIX
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- pheno_vad_MIX[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_MIX_result.RData"))
  
  ### UNK
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- pheno_vad_UNK[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0(trait,"_GeneCentric_Noncoding_STAARO_UNK_result.RData"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  
  
  
  ##### Gene Centric Coding: Burden
  r2_tun_GeneCentric_Coding_Burden  <- rep(0,15)
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_Burden)
  for(k in 1:15){
    prs <- pheno_tuning_Burden[,paste0("Coding_Threshold",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning_Burden)
    r2_tun_GeneCentric_Coding_Burden[k] <- summary(model.prs)$r.square
  }
  
  ### EUR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- pheno_vad_EUR[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_EUR_result.RData"))
  
  ### NonEur
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- pheno_vad_NonEur[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_NonEur_result.RData"))
  
  ### EAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- pheno_vad_EAS[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_EAS_result.RData"))
  
  ### AFR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- pheno_vad_AFR[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_AFR_result.RData"))
  
  ### SAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- pheno_vad_SAS[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_SAS_result.RData"))
  
  ### MIX
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- pheno_vad_MIX[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_MIX_result.RData"))
  
  ### UNK
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- pheno_vad_UNK[,paste0("Coding_Threshold",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  
  save(r2.result_GeneCentric_Coding_Burden, file = paste0(trait,"_GeneCentric_Coding_Burden_UNK_result.RData"))
  
  ##### Gene Centric Noncoding: Burden
  r2_tun_GeneCentric_Noncoding_Burden  <- rep(0,15)
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_Burden)
  for(k in 1:15){
    prs <- pheno_tuning_Burden[,paste0("Noncoding_Threshold",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning_Burden)
    r2_tun_GeneCentric_Noncoding_Burden[k] <- summary(model.prs)$r.square
  }
  
  ### EUR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- pheno_vad_EUR[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_EUR_result.RData"))
  
  ### NonEur
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- pheno_vad_NonEur[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_NonEur_result.RData"))
  
  ### EAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- pheno_vad_EAS[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_EAS_result.RData"))
  
  ### AFR
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- pheno_vad_AFR[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_AFR_result.RData"))
  
  ### SAS
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- pheno_vad_SAS[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_SAS_result.RData"))
  
  ### MIX
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- pheno_vad_MIX[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_MIX_result.RData"))
  
  ### UNK
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- pheno_vad_UNK[,paste0("Noncoding_Threshold",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = paste0(trait,"_GeneCentric_Noncoding_Burden_UNK_result.RData"))
  
}

system("rm All_Tune.txt")
system("rm All_Validation.txt")
system("rm all_phenotypes.RData")