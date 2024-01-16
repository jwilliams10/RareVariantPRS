rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

trait <- "BMI"

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  #Train
  prs_train_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_train.txt"), sep="")
  prs_train_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_prs_train.sscore"))
  prs_train_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_prs_train.sscore"))
  
  prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],prs_train_LDPred2[,-c(1,2,3,4)],prs_train_LASSOSum2[,-c(1,2,3,4)])
  colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
  rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)
  
  ## Merge covariates and y for tuning with the prs_mat
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")
  
  #Tune
  prs_tune_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
  prs_tune_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_prs_tune.sscore"))
  prs_tune_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_prs_tune.sscore"))
  
  prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
  colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
  rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)
  
  ## Merge covariates and y for tuning with the prs_mat
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")
  
  #Validation
  prs_validation_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
  prs_validation_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_prs_validation.sscore"))
  prs_validation_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_prs_validation.sscore"))
  
  prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],prs_validation_LDPred2[,-c(1,2,3,4)],prs_validation_LASSOSum2[,-c(1,2,3,4)])
  colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
  rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)
  
  ## Merge covariates and y for tuning with the prs_mat
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  pheno_validation <- left_join(pheno_validation,prs_validation_all,by = "IID")
  
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  pheno_validation_EUR <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_validation_NonEur <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_validation_UNK <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_validation_SAS <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_validation_MIX <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_validation_AFR <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_validation_EAS <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  
  set.seed(1330)
  
  pheno_tune <- pheno_tune[sample(1:nrow(pheno_tune),size = round(nrow(pheno_validation_NonEur)/2),replace = FALSE),]
  idx_UNK <- sample(1:nrow(pheno_validation_UNK),size = round(nrow(pheno_validation_UNK)/2),replace = FALSE)
  print(summary(idx_UNK))
  idx_SAS <- sample(1:nrow(pheno_validation_SAS),size = round(nrow(pheno_validation_SAS)/2),replace = FALSE)
  idx_MIX <- sample(1:nrow(pheno_validation_MIX),size = round(nrow(pheno_validation_MIX)/2),replace = FALSE)
  idx_AFR <- sample(1:nrow(pheno_validation_AFR),size = round(nrow(pheno_validation_AFR)/2),replace = FALSE)
  idx_EAS <- sample(1:nrow(pheno_validation_EAS),size = round(nrow(pheno_validation_EAS)/2),replace = FALSE)
  
  pheno_tune <- rbind(pheno_tune,pheno_validation_UNK[idx_UNK,],pheno_validation_SAS[idx_SAS,],pheno_validation_MIX[idx_MIX,],pheno_validation_AFR[idx_AFR,],pheno_validation_EAS[idx_EAS,])
  
  pheno_validation_UNK <- pheno_validation_UNK[!((1:nrow(pheno_validation_UNK)) %in% idx_UNK),]
  pheno_validation_SAS <- pheno_validation_SAS[!((1:nrow(pheno_validation_SAS)) %in% idx_SAS),]
  pheno_validation_MIX <- pheno_validation_MIX[!((1:nrow(pheno_validation_MIX)) %in% idx_MIX),]
  pheno_validation_AFR <- pheno_validation_AFR[!((1:nrow(pheno_validation_AFR)) %in% idx_AFR),]
  pheno_validation_EAS <- pheno_validation_EAS[!((1:nrow(pheno_validation_EAS)) %in% idx_EAS),]
  
  pheno_validation_NonEur <- pheno_validation[pheno_validation$IID %in% rbind(pheno_validation_UNK,pheno_validation_SAS,pheno_validation_MIX,pheno_validation_AFR,pheno_validation_EAS)$IID,]
  pheno_validation <- pheno_validation[!(pheno_validation$IID %in% pheno_tune$IID),]
  
  prs_train_all <- pheno_train[!is.na(pheno_train[,trait]),-c(1:27)]
  prs_tune_all <- pheno_tune[!is.na(pheno_tune[,trait]),-c(1:27)]
  prs_validation_all <- pheno_validation[!is.na(pheno_validation[,trait]),-c(1:27)]
  
  ## Null Models
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
  y_tune <- model.null$residual
  
  ##################################
  
  drop <- caret::findLinearCombos(prs_tune_all)$remove
  drop <- names(prs_tune_all)[drop]
  
  prs_train_all = prs_train_all %>% 
    select(-all_of(drop))
  prs_tune_all = prs_tune_all %>% 
    select(-all_of(drop))
  prs_validation_all = prs_validation_all %>% 
    select(-all_of(drop))
  
  mtx <- cor(prs_tune_all)
  drop <- findCorrelation(mtx,cutoff=0.98)
  drop <- names(prs_tune_all)[drop]
  
  prs_train_all = prs_train_all %>% 
    select(-all_of(drop))
  prs_tune_all = prs_tune_all %>% 
    select(-all_of(drop))
  prs_validation_all = prs_validation_all %>% 
    select(-all_of(drop))
  
  
  ## SL
  
  SL.library <- c(
    "SL.glmnet",
    "SL.ridge",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = y_tune, X = prs_tune_all, family = gaussian(),
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library)
  cvsl <- CV.SuperLearner(Y = y_tune, X = prs_tune_all, family = gaussian(),
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library)
  
  best_algorithm <- summary(cvsl)$Table$Algorithm[which.min(summary(cvsl)$Table$Ave)]
  
  a <- predict(sl, prs_validation_all, onlySL = FALSE)
  
  prs_best_validation_sl <- a$pred
  prs_best_validation_glmnet <- a$library.predict[,1]
  prs_best_validation_ridge <- a$library.predict[,2]
  prs_best_validation_glm <- a$library.predict[,3]
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    prs_best_validation <- prs_best_validation_glmnet
  }else if(best_algorithm == "SL.ridge_All"){
    #final
    prs_best_validation <- prs_best_validation_ridge
  }else if(best_algorithm == "SL.glm_All"){
    #final
    prs_best_validation <- prs_best_validation_glm
  }else{
    #final
    prs_best_validation <- prs_best_validation_sl
  }
  
  prs_best_validation <- data.frame(IID = pheno_validation$IID[!is.na(pheno_validation[,trait])],prs = prs_best_validation)
  
  
  a <- predict(sl, prs_train_all, onlySL = FALSE)
  
  prs_best_train_sl <- a$pred
  prs_best_train_glmnet <- a$library.predict[,1]
  prs_best_train_ridge <- a$library.predict[,2]
  prs_best_train_glm <- a$library.predict[,3]
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    prs_best_train <- prs_best_train_glmnet
  }else if(best_algorithm == "SL.ridge_All"){
    #final
    prs_best_train <- prs_best_train_ridge
  }else if(best_algorithm == "SL.glm_All"){
    #final
    prs_best_train <- prs_best_train_glm
  }else{
    #final
    prs_best_train <- prs_best_train_sl
  }
  prs_best_train <- data.frame(IID = pheno_train$IID[!is.na(pheno_train[,trait])],prs = prs_best_train)
  
  
  a <- predict(sl, prs_tune_all, onlySL = FALSE)
  
  prs_best_tune_sl <- a$pred
  prs_best_tune_glmnet <- a$library.predict[,1]
  prs_best_tune_ridge <- a$library.predict[,2]
  prs_best_tune_glm <- a$library.predict[,3]
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    prs_best_tune <- prs_best_tune_glmnet
  }else if(best_algorithm == "SL.ridge_All"){
    #final
    prs_best_tune <- prs_best_tune_ridge
  }else if(best_algorithm == "SL.glm_All"){
    #final
    prs_best_tune <- prs_best_tune_glm
  }else{
    #final
    prs_best_tune <- prs_best_tune_sl
  }
  prs_best_tune <- data.frame(IID = pheno_tune$IID[!is.na(pheno_tune[,trait])],prs = prs_best_tune)
  
  write.table(prs_best_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Train_All_modified.txt"),sep = "\t",row.names = FALSE)
  write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Tune_All_modified.txt"),sep = "\t",row.names = FALSE)
  write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Validation_All_modified.txt"),sep = "\t",row.names = FALSE)
  
  pheno_validation <- left_join(pheno_validation,prs_best_validation)
  
  pheno_validation_EUR <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_validation_NonEur <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_validation_UNK <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_validation_SAS <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_validation_MIX <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_validation_AFR <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_validation_EAS <- pheno_validation[pheno_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_EUR)
  prs <- pheno_validation_EUR[!is.na(pheno_validation_EUR[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_Eur",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur_modified.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_NonEur)
  prs <- pheno_validation_NonEur[!is.na(pheno_validation_NonEur[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_NonEur",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur_modified.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_UNK)
  prs <- pheno_validation_UNK[!is.na(pheno_validation_UNK[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_UNK",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK_modified.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_SAS)
  prs <- pheno_validation_SAS[!is.na(pheno_validation_SAS[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_SAS",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS_modified.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_MIX)
  prs <- pheno_validation_MIX[!is.na(pheno_validation_MIX[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_MIX",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX_modified.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_AFR)
  prs <- pheno_validation_AFR[!is.na(pheno_validation_AFR[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_AFR",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR_modified.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <-  lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation_EAS)
  prs <- pheno_validation_EAS[!is.na(pheno_validation_EAS[,trait]),"prs"]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
  
  ci_result <- boot.ci(boot_r2, type = "perc")
  SL.result <- data.frame(method = "SL_Combined_EAS",
                          r2 = r2,
                          r2_low = ci_result$percent[4],
                          r2_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS_modified.RData"))
  
  
}