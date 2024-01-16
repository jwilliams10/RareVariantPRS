rm(list = ls())
library(data.table)
library(dplyr)
library(boot)

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
  
  prs_mat_train <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_train.txt"), sep="")
  
  prs_mat_tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
  
  prs_mat_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
  
  
  ## Pull in Phenotypes/Covariates 
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")
  
  pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  
  set.seed(1330)
  
  pheno_tuning <- pheno_tuning[sample(1:nrow(pheno_tuning),size = round(nrow(pheno_vad_NonEur)/2),replace = FALSE),]
  idx_UNK <- sample(1:nrow(pheno_vad_UNK),size = round(nrow(pheno_vad_UNK)/2),replace = FALSE)
  idx_SAS <- sample(1:nrow(pheno_vad_SAS),size = round(nrow(pheno_vad_SAS)/2),replace = FALSE)
  idx_MIX <- sample(1:nrow(pheno_vad_MIX),size = round(nrow(pheno_vad_MIX)/2),replace = FALSE)
  idx_AFR <- sample(1:nrow(pheno_vad_AFR),size = round(nrow(pheno_vad_AFR)/2),replace = FALSE)
  idx_EAS <- sample(1:nrow(pheno_vad_EAS),size = round(nrow(pheno_vad_EAS)/2),replace = FALSE)
  
  pheno_tuning <- rbind(pheno_tuning,pheno_vad_UNK[idx_UNK,],pheno_vad_SAS[idx_SAS,],pheno_vad_MIX[idx_MIX,],pheno_vad_AFR[idx_AFR,],pheno_vad_EAS[idx_EAS,])
  
  pheno_vad_UNK <- pheno_vad_UNK[!((1:nrow(pheno_vad_UNK)) %in% idx_UNK),]
  pheno_vad_SAS <- pheno_vad_SAS[!((1:nrow(pheno_vad_SAS)) %in% idx_SAS),]
  pheno_vad_MIX <- pheno_vad_MIX[!((1:nrow(pheno_vad_MIX)) %in% idx_MIX),]
  pheno_vad_AFR <- pheno_vad_AFR[!((1:nrow(pheno_vad_AFR)) %in% idx_AFR),]
  pheno_vad_EAS <- pheno_vad_EAS[!((1:nrow(pheno_vad_EAS)) %in% idx_EAS),]
  
  pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% rbind(pheno_vad_UNK,pheno_vad_SAS,pheno_vad_MIX,pheno_vad_AFR,pheno_vad_EAS)$IID,]
  
  
  #calculate R2 for each of the tuning dataset
  # This is done by regressing the residuals of the model with all covariates against the prs
  r2_tun_vec <- rep(0,length(pthres))
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning)
  for(k in 1:length(pthres)){
    prs <- pheno_tuning[!is.na(pheno_tuning[,trait]),paste0("p_value_",k)]
    model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
    r2_tun_vec[k] <- summary(model.prs)$r.square
  }
  #find best p-value threshold
  idx <- which.max(r2_tun_vec)
  #write down best prs in the prs folder/ save it and store it
  prs_train_max <- pheno_train[,c("IID","FID",paste0("p_value_",idx))]
  colnames(prs_train_max) <- c("IID","FID","prs")
  write.table(prs_train_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_train_best_modified.txt"),row.names = F)
  
  prs_tune_max <- pheno_tuning[,c("IID","FID",paste0("p_value_",idx))]
  colnames(prs_tune_max) <- c("IID","FID","prs")
  write.table(prs_tune_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_tune_best_modified.txt"),row.names = F)
  
  prs_vad_max <- pheno_vad[,c("IID","FID",paste0("p_value_",idx))]
  colnames(prs_vad_max) <- c("IID","FID","prs")
  write.table(prs_vad_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_validation_best_modified.txt"),row.names = F)
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_EUR",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EUR_modified.RData")) 
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_NonEur",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_NonEur_modified.RData")) 
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_UNK",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_UNK_modified.RData")) 
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_SAS",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_SAS_modified.RData")) 
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_MIX",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_MIX_modified.RData"))
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_AFR",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_AFR_modified.RData"))
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),paste0("p_value_",idx)]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  r2.result <- data.frame(method = "CT_EAS",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  ## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
  ct.result <- list(r2.result,r2_tun_vec)
  save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EAS_modified.RData"))
  
  
}