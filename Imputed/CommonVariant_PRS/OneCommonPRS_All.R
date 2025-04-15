rm(list = ls())
library(data.table)
library(dplyr)
library(RISCA)
library(boot)
library(stringr)
library(caret)
library(ranger)
library(glmnet)

time <- system.time({
  continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
  
  binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")
  
  trait <- as.numeric(commandArgs(TRUE)[1])
  
  trait <- c(continuous_traits,binary_traits)[trait]
  
  CT_prs_all_train <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_train.txt"), sep="")
  LDpred2_prs_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train.sscore"), header=FALSE, comment.char="#")
  LASSOsum2_prs_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train.sscore"), header=FALSE, comment.char="#")
  
  prs_train_all <- cbind(LDpred2_prs_train[,1:2],CT_prs_all_train[,-c(1,2)],LDpred2_prs_train[,-c(1,2,3,4)],LASSOsum2_prs_train[,-c(1,2,3,4)])
  colnames(prs_train_all) <- c("FID","IID",paste0("CT_",colnames(CT_prs_all_train[,-c(1,2)]),"_",trait),paste0("LDPred2_",colnames(LDpred2_prs_train[,-c(1,2,3,4)]),"_",trait),paste0("LASSOSum2_",colnames(LASSOsum2_prs_train[,-c(1,2,3,4)]),"_",trait))
  rm(CT_prs_all_train);rm(LDpred2_prs_train);rm(LASSOsum2_prs_train)
  
  CT_prs_all_tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
  LDpred2_prs_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune.sscore"), header=FALSE, comment.char="#")
  LASSOsum2_prs_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune.sscore"), header=FALSE, comment.char="#")
  
  prs_tune_all <- cbind(LDpred2_prs_tune[,1:2],CT_prs_all_tune[,-c(1,2)],LDpred2_prs_tune[,-c(1,2,3,4)],LASSOsum2_prs_tune[,-c(1,2,3,4)])
  colnames(prs_tune_all) <- c("FID","IID",paste0("CT_",colnames(CT_prs_all_tune[,-c(1,2)]),"_",trait),paste0("LDPred2_",colnames(LDpred2_prs_tune[,-c(1,2,3,4)]),"_",trait),paste0("LASSOSum2_",colnames(LASSOsum2_prs_tune[,-c(1,2,3,4)]),"_",trait))
  rm(CT_prs_all_tune);rm(LDpred2_prs_tune);rm(LASSOsum2_prs_tune)
  
  CT_prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
  LDpred2_prs_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation.sscore"), header=FALSE, comment.char="#")
  LASSOsum2_prs_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation.sscore"), header=FALSE, comment.char="#")
  
  prs_validation_all <- cbind(LDpred2_prs_validation[,1:2],CT_prs_all_validation[,-c(1,2)],LDpred2_prs_validation[,-c(1,2,3,4)],LASSOsum2_prs_validation[,-c(1,2,3,4)])
  colnames(prs_validation_all) <- c("FID","IID",paste0("CT_",colnames(CT_prs_all_validation[,-c(1,2)]),"_",trait),paste0("LDPred2_",colnames(LDpred2_prs_validation[,-c(1,2,3,4)]),"_",trait),paste0("LASSOSum2_",colnames(LASSOsum2_prs_validation[,-c(1,2,3,4)]),"_",trait))
  rm(CT_prs_all_validation);rm(LDpred2_prs_validation);rm(LASSOsum2_prs_validation)
  
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  pheno_train <- left_join(pheno_train,prs_train_all)
  prs_train_all <- pheno_train[,str_detect(colnames(pheno_train),paste0("_",trait))]
  
  ## Merge covariates and y for tuning with the prs_mat
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  pheno_tune <- left_join(pheno_tune,prs_tune_all)
  prs_tune_all <- pheno_tune[,str_detect(colnames(pheno_tune),paste0("_",trait))]
  
  ## Merge covariates and y for tuning with the prs_mat
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  pheno_validation <- left_join(pheno_validation,prs_validation_all)
  prs_validation_all <- pheno_validation[,str_detect(colnames(pheno_validation),paste0("_",trait))]
  
  
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
  
  
  Ensemble_Function_Continuous <- function(x,y){
    x <- as.matrix(x[!is.na(y),])
    y <- y[!is.na(y)]
    
    lasso_mod <- cv.glmnet(x,y,family = "gaussian",alpha = 1,type.measure = "mse",nfold = 10)
    ridge_mod <- cv.glmnet(x,y,family = "gaussian",alpha = 0,type.measure = "mse",nfold = 10)
    
    lasso_prediction_x <- predict(lasso_mod, x)
    ridge_prediction_x <- predict(ridge_mod, x)
    
    ensemble_mod <- lm(y~.,data = data.frame(lasso_prediction_x,ridge_prediction_x))
    
    ensemble_prediction_x <- ensemble_mod$fitted
    
    coefficients_x <- coef(lm(y~.,data.frame(y = ensemble_prediction_x,x)))
    return(list(Coefficients = coefficients_x))
  }
  Ensemble_Function_Binary<- function(x,y){
    x <- as.matrix(x[!is.na(y),])
    y <- y[!is.na(y)]
    
    lasso_mod <- cv.glmnet(x,y,family = "binomial",alpha = 1,type.measure = "auc")
    ridge_mod <- cv.glmnet(x,y,family = "binomial",alpha = 0,type.measure = "auc")
    
    lasso_prediction_x <- predict(lasso_mod, x,type = "link")
    ridge_prediction_x <- predict(ridge_mod, x,type = "link")
    
    ensemble_mod <- glm(y~.,data = data.frame(lasso_prediction_x,ridge_prediction_x),family = binomial())
    ensemble_prediction_x <- predict(ensemble_mod,data.frame(lasso_prediction_x,ridge_prediction_x),type = "link")
    
    coefficients_x <- coef(lm(y~.,data.frame(y = ensemble_prediction_x,x)))
    return(list(Coefficients = coefficients_x))
  }
  Ensemble_Function <- function(x,y,family = c("continuous","binary")){
    if(family == "continuous"){
      return(Ensemble_Function_Continuous(x,y))
    }else{
      return(Ensemble_Function_Binary(x,y))
    }
  }
  
  if(trait %in% binary_traits){
    
    Results <- Ensemble_Function(x = prs_tune_all,y = pheno_tune[,trait],family = "binary")
    Results$Coefficients[is.na(Results$Coefficients)] <- 0
    save(Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_Coefficients.csv"))
    PRS_Train <- as.matrix(pheno_train[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    PRS_Tune <- as.matrix(pheno_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    PRS_Validation <- as.matrix(pheno_validation[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    
    PRS_Train <- data.frame(IID = pheno_train$IID,PRS = PRS_Train)
    write.csv(PRS_Train,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Train.csv"),row.names = FALSE)
    PRS_Tune <- data.frame(IID = pheno_tune$IID,PRS = PRS_Tune)
    write.csv(PRS_Tune,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"),row.names = FALSE)
    PRS_Validation <- data.frame(IID = pheno_validation$IID,PRS = PRS_Validation)
    write.csv(PRS_Validation,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"),row.names = FALSE)
    
    if(trait %in% c("Breast","Prostate")){
      fill <- "withoutsex"
      confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }else{
      fill <- "withsex"
      confounders <- paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_validation_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),trait)],PRS_Validation)
    pheno_validation_adjusted <- pheno_validation_raw
    tmp <- data.frame(y = pheno_validation_adjusted[,"PRS"],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      pheno_validation_adjusted[,"PRS"] <- 0
    }else{
      pheno_validation_adjusted[,"PRS"] <- R/sqrt(y_hat)
    }
    
    pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    pheno_validation_raw_EUR[,"PRS"] <- scale(pheno_validation_raw_EUR[,"PRS"])
    pheno_validation_raw_SAS[,"PRS"] <- scale(pheno_validation_raw_SAS[,"PRS"])
    pheno_validation_raw_AMR[,"PRS"] <- scale(pheno_validation_raw_AMR[,"PRS"])
    pheno_validation_raw_AFR[,"PRS"] <- scale(pheno_validation_raw_AFR[,"PRS"])
    pheno_validation_raw_EAS[,"PRS"] <- scale(pheno_validation_raw_EAS[,"PRS"])
    
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
      return(c(result))
    }
    
    AUC_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 10000)
    beta_raw_EUR_boot <- boot_beta$t
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    
    AUC_validation_raw_EUR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_raw_EUR[!is.na(pheno_validation_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_EUR, statistic = AUC_Boot, R = 10000)
    AUC_raw_EUR_boot <- boot_AUC$t
    AUC_se_validation_raw_EUR <- sd(boot_AUC$t)
    
    beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 10000)
    beta_raw_SAS_boot <- boot_beta$t
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    
    AUC_validation_raw_SAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_raw_SAS[!is.na(pheno_validation_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_SAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_SAS_boot <- boot_AUC$t
    AUC_se_validation_raw_SAS <- sd(boot_AUC$t)
    
    beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 10000)
    beta_raw_AMR_boot <- boot_beta$t
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    
    AUC_validation_raw_AMR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_raw_AMR[!is.na(pheno_validation_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_AMR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AMR_boot <- boot_AUC$t
    AUC_se_validation_raw_AMR <- sd(boot_AUC$t)
    
    beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 10000)
    beta_raw_AFR_boot <- boot_beta$t
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    
    AUC_validation_raw_AFR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_raw_AFR[!is.na(pheno_validation_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_AFR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AFR_boot <- boot_AUC$t
    AUC_se_validation_raw_AFR <- sd(boot_AUC$t)
    
    beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 10000)
    beta_raw_EAS_boot <- boot_beta$t
    beta_se_validation_raw_EAS <- sd(boot_beta$t)
    
    AUC_validation_raw_EAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_raw_EAS[!is.na(pheno_validation_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_EAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_EAS_boot <- boot_AUC$t
    AUC_se_validation_raw_EAS <- sd(boot_AUC$t)
    
    beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EUR_boot <- boot_beta$t
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_EUR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_adjusted_EUR[!is.na(pheno_validation_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_EUR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EUR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_EUR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_SAS_boot <- boot_beta$t
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    AUC_validation_adjusted_SAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_adjusted_SAS[!is.na(pheno_validation_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_SAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_SAS_boot <- boot_AUC$t
    AUC_se_validation_adjusted_SAS <- sd(boot_AUC$t)
    
    beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AMR_boot <- boot_beta$t
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_AMR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_adjusted_AMR[!is.na(pheno_validation_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_AMR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AMR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_AMR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AFR_boot <- boot_beta$t
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_AFR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_adjusted_AFR[!is.na(pheno_validation_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_AFR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AFR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_AFR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EAS_boot <- boot_beta$t
    beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    AUC_validation_adjusted_EAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = pheno_validation_adjusted_EAS[!is.na(pheno_validation_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_EAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EAS_boot <- boot_AUC$t
    AUC_se_validation_adjusted_EAS <- sd(boot_AUC$t)
    
    RICE_CV <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                 beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                 beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                                 AUC_raw = c(AUC_validation_raw_EUR,AUC_validation_raw_SAS,AUC_validation_raw_AMR,AUC_validation_raw_AFR,AUC_validation_raw_EAS),
                                 AUC_se_raw = c(AUC_se_validation_raw_EUR,AUC_se_validation_raw_SAS,AUC_se_validation_raw_AMR,AUC_se_validation_raw_AFR,AUC_se_validation_raw_EAS),
                                 beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                 beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                                 AUC_adjusted = c(AUC_validation_adjusted_EUR,AUC_validation_adjusted_SAS,AUC_validation_adjusted_AMR,AUC_validation_adjusted_AFR,AUC_validation_adjusted_EAS),
                                 AUC_se_adjusted = c(AUC_se_validation_adjusted_EUR,AUC_se_validation_adjusted_SAS,AUC_se_validation_adjusted_AMR,AUC_se_validation_adjusted_AFR,AUC_se_validation_adjusted_EAS))
    
    RICE_CV_Boot <- data.frame(trait = trait,beta_raw_EUR_boot,AUC_raw_EUR_boot,beta_raw_SAS_boot,AUC_raw_SAS_boot,
                                  beta_raw_AMR_boot,AUC_raw_AMR_boot,beta_raw_AFR_boot,AUC_raw_AFR_boot,
                                  beta_raw_EAS_boot,AUC_raw_EAS_boot,beta_adjusted_EUR_boot,AUC_adjusted_EUR_boot,
                                  beta_adjusted_SAS_boot,AUC_adjusted_SAS_boot,beta_adjusted_AMR_boot,AUC_adjusted_AMR_boot,
                                  beta_adjusted_AFR_boot,AUC_adjusted_AFR_boot,beta_adjusted_EAS_boot,AUC_adjusted_EAS_boot)
  }else{
    ## Null Models
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_train)
    y_train <- model.null$residual
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
    y_tune <- model.null$residual
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
    y_validation <- model.null$residual
    
    pheno_tune$y_tune <- NA
    pheno_tune$y_tune[!is.na(pheno_tune[,trait])] <- y_tune
    
    pheno_validation$y_validation <- NA
    pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- y_validation
    
    Results <- Ensemble_Function(x = prs_tune_all,y = pheno_tune[,"y_tune"],family = "continuous")
    Results$Coefficients[is.na(Results$Coefficients)] <- 0
    save(Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_Coefficients.csv"))
    PRS_Train <- as.matrix(pheno_train[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    PRS_Tune <- as.matrix(pheno_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    PRS_Validation <- as.matrix(pheno_validation[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    
    PRS_Train <- data.frame(IID = pheno_train$IID,PRS = PRS_Train)
    write.csv(PRS_Train,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Train.csv"),row.names = FALSE)
    PRS_Tune <- data.frame(IID = pheno_tune$IID,PRS = PRS_Tune)
    write.csv(PRS_Tune,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"),row.names = FALSE)
    PRS_Validation <- data.frame(IID = pheno_validation$IID,PRS = PRS_Validation)
    write.csv(PRS_Validation,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"),row.names = FALSE)
    
    pheno_validation_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),"y_validation")],PRS_Validation)
    pheno_validation_adjusted <- pheno_validation_raw
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_validation_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),"y_validation")],PRS_Validation)
    pheno_validation_adjusted <- pheno_validation_raw
    
    tmp <- data.frame(y = pheno_validation_adjusted[,"PRS"],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(sqrt(y_hat)) == 0){
      pheno_validation_adjusted[,"PRS"] <- 0
    }else{
      pheno_validation_adjusted[,"PRS"] <- R/sqrt(y_hat)
    }
    
    pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
    pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
    pheno_validation_raw_AMR$y_validation <- scale(pheno_validation_raw_AMR$y_validation)
    pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)
    pheno_validation_raw_EAS$y_validation <- scale(pheno_validation_raw_EAS$y_validation)
    
    pheno_validation_adjusted_EUR$y_validation <- scale(pheno_validation_adjusted_EUR$y_validation)
    pheno_validation_adjusted_SAS$y_validation <- scale(pheno_validation_adjusted_SAS$y_validation)
    pheno_validation_adjusted_AMR$y_validation <- scale(pheno_validation_adjusted_AMR$y_validation)
    pheno_validation_adjusted_AFR$y_validation <- scale(pheno_validation_adjusted_AFR$y_validation)
    pheno_validation_adjusted_EAS$y_validation <- scale(pheno_validation_adjusted_EAS$y_validation)
    
    pheno_validation_raw_EUR[,"PRS"] <- scale(pheno_validation_raw_EUR[,"PRS"])
    pheno_validation_raw_SAS[,"PRS"] <- scale(pheno_validation_raw_SAS[,"PRS"])
    pheno_validation_raw_AMR[,"PRS"] <- scale(pheno_validation_raw_AMR[,"PRS"])
    pheno_validation_raw_AFR[,"PRS"] <- scale(pheno_validation_raw_AFR[,"PRS"])
    pheno_validation_raw_EAS[,"PRS"] <- scale(pheno_validation_raw_EAS[,"PRS"])
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = boot_data))[2]
      return(c(result))
    }
    
    R2_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = boot_data))$r.squared
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_EUR))[2]
    boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 10000)
    beta_raw_EUR_boot <- boot_beta$t
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    
    R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_EUR))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 10000)
    R2_raw_EUR_boot <- boot_R2$t
    R2_se_validation_raw_EUR <- sd(boot_R2$t)
    
    beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_SAS))[2]
    boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 10000)
    beta_raw_SAS_boot <- boot_beta$t
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    
    R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_SAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 10000)
    R2_raw_SAS_boot <- boot_R2$t
    R2_se_validation_raw_SAS <- sd(boot_R2$t)
    
    beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AMR))[2]
    boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 10000)
    beta_raw_AMR_boot <- boot_beta$t
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    
    R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AMR))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 10000)
    R2_raw_AMR_boot <- boot_R2$t
    R2_se_validation_raw_AMR <- sd(boot_R2$t)
    
    beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AFR))[2]
    boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 10000)
    beta_raw_AFR_boot <- boot_beta$t
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    
    R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AFR))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 10000)
    R2_raw_AFR_boot <- boot_R2$t
    R2_se_validation_raw_AFR <- sd(boot_R2$t)
    
    beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_EAS))[2]
    boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 10000)
    beta_raw_EAS_boot <- boot_beta$t
    beta_se_validation_raw_EAS <- sd(boot_beta$t)
    
    R2_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_EAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_EAS, statistic = R2_Boot, R = 10000)
    R2_raw_EAS_boot <- boot_R2$t
    R2_se_validation_raw_EAS <- sd(boot_R2$t)
    
    beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_EUR))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EUR_boot <- boot_beta$t
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_EUR))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 10000)
    R2_adjusted_EUR_boot <- boot_R2$t
    R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
    
    beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_SAS))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_SAS_boot <- boot_beta$t
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_SAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_SAS_boot <- boot_R2$t
    R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
    
    beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AMR))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AMR_boot <- boot_beta$t
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AMR))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AMR_boot <- boot_R2$t
    R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
    
    beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AFR))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AFR_boot <- boot_beta$t
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AFR))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AFR_boot <- boot_R2$t
    R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
    
    beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_EAS))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EAS_boot <- boot_beta$t
    beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_EAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_EAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_EAS_boot <- boot_R2$t
    R2_se_validation_adjusted_EAS <- sd(boot_R2$t)
    
    RICE_CV <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                 beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                 beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                                 R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                                 R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                                 beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                 beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                                 R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                                 R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))
    
    RICE_CV_Boot <- data.frame(trait = trait,beta_raw_EUR_boot,R2_raw_EUR_boot,beta_raw_SAS_boot,R2_raw_SAS_boot,
                                  beta_raw_AMR_boot,R2_raw_AMR_boot,beta_raw_AFR_boot,R2_raw_AFR_boot,
                                  beta_raw_EAS_boot,R2_raw_EAS_boot,beta_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                                  beta_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                                  beta_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_adjusted_EAS_boot,R2_adjusted_EAS_boot)
  }
  write.csv(RICE_CV,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"Best_Betas.csv"),row.names = FALSE)
  write.csv(RICE_CV_Boot,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_Bootstraps.csv"),row.names = FALSE)
})[3]
save(time, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_Time.RData"))