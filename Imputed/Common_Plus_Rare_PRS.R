rm(list = ls())
library(readr)
library(dplyr)
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(stringr)
library(glmnet)
library(RISCA)

time <- system.time({
  
  continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
  binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")
  trait <- c(continuous_traits,binary_traits)[as.numeric(commandArgs(TRUE)[1])]
  
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  CV_PRS_Tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"))
  colnames(CV_PRS_Tune) <- c("IID","CV_PRS")
  pheno_tune <- inner_join(pheno_tune,CV_PRS_Tune)
  RV_PRS_Tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Tune.csv"))
  colnames(RV_PRS_Tune) <- c("IID","RV_PRS")
  pheno_tune <- inner_join(pheno_tune,RV_PRS_Tune)
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  if(trait %in% continuous_traits){
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
    pheno_tune$y_tune <- NA
    pheno_tune$y_tune[!is.na(pheno_tune[,trait])] <- model.null$residual
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
    pheno_validation$y_validation <- NA
    pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
    
    RICE_Model <- lm(y_tune ~ CV_PRS + RV_PRS,data = pheno_tune)
    pheno_validation$PRS <- predict(RICE_Model,pheno_validation)
    
    CV_RV_PRS_raw <- pheno_validation
    CV_RV_PRS_adjusted <- pheno_validation
    
    for(i in c("CV_PRS","RV_PRS","PRS")){
      tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
      mod <- lm(y~.,data = tmp)
      R <- mod$residuals
      tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
      mod <- lm(y~.,data = tmp)
      y_hat <- predict(mod,tmp)
      if(sum(y_hat < 0) > 0){
        mod <- lm(y~1,data = tmp)
        y_hat <- predict(mod,tmp)
      }
      if(sum(sqrt(y_hat)) == 0){
        CV_RV_PRS_adjusted[,i] <- 0
      }else{
        CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
      }
    }
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    CV_RV_PRS_raw_EUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    CV_RV_PRS_raw_SAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    CV_RV_PRS_raw_AMR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    CV_RV_PRS_raw_AFR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    CV_RV_PRS_raw_EAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    CV_RV_PRS_adjusted_EAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    CV_RV_PRS_raw_EUR$y_validation <- scale(CV_RV_PRS_raw_EUR$y_validation)
    CV_RV_PRS_raw_SAS$y_validation <- scale(CV_RV_PRS_raw_SAS$y_validation)
    CV_RV_PRS_raw_AMR$y_validation <- scale(CV_RV_PRS_raw_AMR$y_validation)
    CV_RV_PRS_raw_AFR$y_validation <- scale(CV_RV_PRS_raw_AFR$y_validation)
    CV_RV_PRS_raw_EAS$y_validation <- scale(CV_RV_PRS_raw_EAS$y_validation)
    
    CV_RV_PRS_adjusted_EUR$y_validation <- scale(CV_RV_PRS_adjusted_EUR$y_validation)
    CV_RV_PRS_adjusted_SAS$y_validation <- scale(CV_RV_PRS_adjusted_SAS$y_validation)
    CV_RV_PRS_adjusted_AMR$y_validation <- scale(CV_RV_PRS_adjusted_AMR$y_validation)
    CV_RV_PRS_adjusted_AFR$y_validation <- scale(CV_RV_PRS_adjusted_AFR$y_validation)
    CV_RV_PRS_adjusted_EAS$y_validation <- scale(CV_RV_PRS_adjusted_EAS$y_validation)
    
    CV_RV_PRS_raw_EUR$RV_PRS <- scale(CV_RV_PRS_raw_EUR$RV_PRS)
    CV_RV_PRS_raw_SAS$RV_PRS <- scale(CV_RV_PRS_raw_SAS$RV_PRS)
    CV_RV_PRS_raw_AMR$RV_PRS <- scale(CV_RV_PRS_raw_AMR$RV_PRS)
    CV_RV_PRS_raw_AFR$RV_PRS <- scale(CV_RV_PRS_raw_AFR$RV_PRS)
    CV_RV_PRS_raw_EAS$RV_PRS <- scale(CV_RV_PRS_raw_EAS$RV_PRS)
    
    CV_RV_PRS_raw_EUR$CV_PRS <- scale(CV_RV_PRS_raw_EUR$CV_PRS)
    CV_RV_PRS_raw_SAS$CV_PRS <- scale(CV_RV_PRS_raw_SAS$CV_PRS)
    CV_RV_PRS_raw_AMR$CV_PRS <- scale(CV_RV_PRS_raw_AMR$CV_PRS)
    CV_RV_PRS_raw_AFR$CV_PRS <- scale(CV_RV_PRS_raw_AFR$CV_PRS)
    CV_RV_PRS_raw_EAS$CV_PRS <- scale(CV_RV_PRS_raw_EAS$CV_PRS)
    
    CV_RV_PRS_raw_EUR$PRS <- scale(CV_RV_PRS_raw_EUR$PRS)
    CV_RV_PRS_raw_SAS$PRS <- scale(CV_RV_PRS_raw_SAS$PRS)
    CV_RV_PRS_raw_AMR$PRS <- scale(CV_RV_PRS_raw_AMR$PRS)
    CV_RV_PRS_raw_AFR$PRS <- scale(CV_RV_PRS_raw_AFR$PRS)
    CV_RV_PRS_raw_EAS$PRS <- scale(CV_RV_PRS_raw_EAS$PRS)
    
    
    Beta_CV_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(lm(y_validation~CV_PRS + RV_PRS,data = boot_data))[2]
      return(c(result))
    }
    
    Beta_RV_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(lm(y_validation~CV_PRS + RV_PRS,data = boot_data))[3]
      return(c(result))
    }
    
    R2_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- summary(lm(y_validation~PRS,data = boot_data))$r.squared
      return(c(result))
    }
    
    beta_CV_validation_raw_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_EUR))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_EUR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_EUR_boot <- boot_beta$t
    beta_CV_se_validation_raw_EUR <- sd(boot_beta$t)
    
    beta_RV_validation_raw_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_EUR))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_EUR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_EUR_boot <- boot_beta$t
    beta_RV_se_validation_raw_EUR <- sd(boot_beta$t)
    
    R2_validation_raw_EUR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_EUR))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_raw_EUR, statistic = R2_Boot, R = 10000)
    R2_raw_EUR_boot <- boot_R2$t
    R2_se_validation_raw_EUR <- sd(boot_R2$t)
    
    beta_CV_validation_raw_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_SAS))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_SAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_SAS_boot <- boot_beta$t
    beta_CV_se_validation_raw_SAS <- sd(boot_beta$t)
    
    beta_RV_validation_raw_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_SAS))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_SAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_SAS_boot <- boot_beta$t
    beta_RV_se_validation_raw_SAS <- sd(boot_beta$t)
    
    R2_validation_raw_SAS <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_SAS))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_raw_SAS, statistic = R2_Boot, R = 10000)
    R2_raw_SAS_boot <- boot_R2$t
    R2_se_validation_raw_SAS <- sd(boot_R2$t)
    
    beta_CV_validation_raw_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AMR))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_AMR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_AMR_boot <- boot_beta$t
    beta_CV_se_validation_raw_AMR <- sd(boot_beta$t)
    
    beta_RV_validation_raw_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AMR))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_AMR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_AMR_boot <- boot_beta$t
    beta_RV_se_validation_raw_AMR <- sd(boot_beta$t)
    
    R2_validation_raw_AMR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_AMR))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_raw_AMR, statistic = R2_Boot, R = 10000)
    R2_raw_AMR_boot <- boot_R2$t
    R2_se_validation_raw_AMR <- sd(boot_R2$t)
    
    beta_CV_validation_raw_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AFR))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_AFR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_AFR_boot <- boot_beta$t
    beta_CV_se_validation_raw_AFR <- sd(boot_beta$t)
    
    beta_RV_validation_raw_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AFR))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_AFR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_AFR_boot <- boot_beta$t
    beta_RV_se_validation_raw_AFR <- sd(boot_beta$t)
    
    R2_validation_raw_AFR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_AFR))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_raw_AFR, statistic = R2_Boot, R = 10000)
    R2_raw_AFR_boot <- boot_R2$t
    R2_se_validation_raw_AFR <- sd(boot_R2$t)
    
    beta_CV_validation_raw_EAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_EAS))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_EAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_EAS_boot <- boot_beta$t
    beta_CV_se_validation_raw_EAS <- sd(boot_beta$t)
    
    beta_RV_validation_raw_EAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_EAS))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_EAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_EAS_boot <- boot_beta$t
    beta_RV_se_validation_raw_EAS <- sd(boot_beta$t)
    
    R2_validation_raw_EAS <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_EAS))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_raw_EAS, statistic = R2_Boot, R = 10000)
    R2_raw_EAS_boot <- boot_R2$t
    R2_se_validation_raw_EAS <- sd(boot_R2$t)
    
    beta_CV_validation_adjusted_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_EUR_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_EUR_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    R2_validation_adjusted_EUR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_EUR))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = R2_Boot, R = 10000)
    R2_adjusted_EUR_boot <- boot_R2$t
    R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
    
    beta_CV_validation_adjusted_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_SAS_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_SAS_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_SAS <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_SAS))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_SAS_boot <- boot_R2$t
    R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
    
    beta_CV_validation_adjusted_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AMR))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_AMR_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AMR))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_AMR_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AMR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_AMR))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AMR_boot <- boot_R2$t
    R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
    
    beta_CV_validation_adjusted_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_AFR_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_AFR_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AFR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_AFR))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AFR_boot <- boot_R2$t
    R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
    
    beta_CV_validation_adjusted_EAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_EAS))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_EAS_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_EAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_EAS))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_EAS_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_EAS <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_EAS))$r.squared
    boot_R2 <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_EAS_boot <- boot_R2$t
    R2_se_validation_adjusted_EAS <- sd(boot_R2$t)
    
    CV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                 beta_raw = c(beta_CV_validation_raw_EUR,beta_CV_validation_raw_SAS,beta_CV_validation_raw_AMR,beta_CV_validation_raw_AFR,beta_CV_validation_raw_EAS), 
                                 beta_se_raw = c(beta_CV_se_validation_raw_EUR,beta_CV_se_validation_raw_SAS,beta_CV_se_validation_raw_AMR,beta_CV_se_validation_raw_AFR,beta_CV_se_validation_raw_EAS), 
                                 R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                                 R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                                 beta_adjusted = c(beta_CV_validation_adjusted_EUR,beta_CV_validation_adjusted_SAS,beta_CV_validation_adjusted_AMR,beta_CV_validation_adjusted_AFR,beta_CV_validation_adjusted_EAS), 
                                 beta_se_adjusted = c(beta_CV_se_validation_adjusted_EUR,beta_CV_se_validation_adjusted_SAS,beta_CV_se_validation_adjusted_AMR,beta_CV_se_validation_adjusted_AFR,beta_CV_se_validation_adjusted_EAS), 
                                 R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                                 R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))
    
    CV_Boot_Results <- data.frame(trait = trait,beta_CV_raw_EUR_boot,R2_raw_EUR_boot,beta_CV_raw_SAS_boot,R2_raw_SAS_boot,
                                  beta_CV_raw_AMR_boot,R2_raw_AMR_boot,beta_CV_raw_AFR_boot,R2_raw_AFR_boot,
                                  beta_CV_raw_EAS_boot,R2_raw_EAS_boot,beta_CV_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                                  beta_CV_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_CV_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                                  beta_CV_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_CV_adjusted_EAS_boot,R2_adjusted_EAS_boot)
    
    RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                 beta_raw = c(beta_RV_validation_raw_EUR,beta_RV_validation_raw_SAS,beta_RV_validation_raw_AMR,beta_RV_validation_raw_AFR,beta_RV_validation_raw_EAS), 
                                 beta_se_raw = c(beta_RV_se_validation_raw_EUR,beta_RV_se_validation_raw_SAS,beta_RV_se_validation_raw_AMR,beta_RV_se_validation_raw_AFR,beta_RV_se_validation_raw_EAS), 
                                 R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                                 R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                                 beta_adjusted = c(beta_RV_validation_adjusted_EUR,beta_RV_validation_adjusted_SAS,beta_RV_validation_adjusted_AMR,beta_RV_validation_adjusted_AFR,beta_RV_validation_adjusted_EAS), 
                                 beta_se_adjusted = c(beta_RV_se_validation_adjusted_EUR,beta_RV_se_validation_adjusted_SAS,beta_RV_se_validation_adjusted_AMR,beta_RV_se_validation_adjusted_AFR,beta_RV_se_validation_adjusted_EAS), 
                                 R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                                 R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))
    
    RV_Boot_Results <- data.frame(trait = trait,beta_RV_raw_EUR_boot,R2_raw_EUR_boot,beta_RV_raw_SAS_boot,R2_raw_SAS_boot,
                                  beta_RV_raw_AMR_boot,R2_raw_AMR_boot,beta_RV_raw_AFR_boot,R2_raw_AFR_boot,
                                  beta_RV_raw_EAS_boot,R2_raw_EAS_boot,beta_RV_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                                  beta_RV_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_RV_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                                  beta_RV_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_RV_adjusted_EAS_boot,R2_adjusted_EAS_boot)
    
  }else{
    RICE_Model <- glm(as.formula(paste0(trait,"~ CV_PRS + RV_PRS")),data = pheno_tune,family = binomial())
    pheno_validation$PRS <- predict(RICE_Model,pheno_validation,type = "link")
    
    CV_RV_PRS_raw <- pheno_validation
    CV_RV_PRS_adjusted <- pheno_validation
    
    for(i in c("CV_PRS","RV_PRS","PRS")){
      tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
      mod <- lm(y~.,data = tmp)
      R <- mod$residuals
      tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
      mod <- lm(y~.,data = tmp)
      y_hat <- predict(mod,tmp)
      if(sum(y_hat < 0) > 0){
        mod <- lm(y~1,data = tmp)
        y_hat <- predict(mod,tmp)
      }
      if(sum(sqrt(y_hat)) == 0){
        CV_RV_PRS_adjusted[,i] <- 0
      }else{
        CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
      }
    }
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    CV_RV_PRS_raw_EUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    CV_RV_PRS_raw_SAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    CV_RV_PRS_raw_AMR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    CV_RV_PRS_raw_AFR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    CV_RV_PRS_raw_EAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    CV_RV_PRS_adjusted_EAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    CV_RV_PRS_raw_EUR$RV_PRS <- scale(CV_RV_PRS_raw_EUR$RV_PRS)
    CV_RV_PRS_raw_SAS$RV_PRS <- scale(CV_RV_PRS_raw_SAS$RV_PRS)
    CV_RV_PRS_raw_AMR$RV_PRS <- scale(CV_RV_PRS_raw_AMR$RV_PRS)
    CV_RV_PRS_raw_AFR$RV_PRS <- scale(CV_RV_PRS_raw_AFR$RV_PRS)
    CV_RV_PRS_raw_EAS$RV_PRS <- scale(CV_RV_PRS_raw_EAS$RV_PRS)
    
    CV_RV_PRS_raw_EUR$CV_PRS <- scale(CV_RV_PRS_raw_EUR$CV_PRS)
    CV_RV_PRS_raw_SAS$CV_PRS <- scale(CV_RV_PRS_raw_SAS$CV_PRS)
    CV_RV_PRS_raw_AMR$CV_PRS <- scale(CV_RV_PRS_raw_AMR$CV_PRS)
    CV_RV_PRS_raw_AFR$CV_PRS <- scale(CV_RV_PRS_raw_AFR$CV_PRS)
    CV_RV_PRS_raw_EAS$CV_PRS <- scale(CV_RV_PRS_raw_EAS$CV_PRS)
    
    CV_RV_PRS_raw_EUR$PRS <- scale(CV_RV_PRS_raw_EUR$PRS)
    CV_RV_PRS_raw_SAS$PRS <- scale(CV_RV_PRS_raw_SAS$PRS)
    CV_RV_PRS_raw_AMR$PRS <- scale(CV_RV_PRS_raw_AMR$PRS)
    CV_RV_PRS_raw_AFR$PRS <- scale(CV_RV_PRS_raw_AFR$PRS)
    CV_RV_PRS_raw_EAS$PRS <- scale(CV_RV_PRS_raw_EAS$PRS)
    
    
    if(trait %in% c("Breast","Prostate")){
      confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }else{
      confounders <-  paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }
    
    Beta_CV_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
      return(c(result))
    }
    
    Beta_RV_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[3]
      return(c(result))
    }
    
    AUC_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
      return(c(result))
    }
    
    beta_CV_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_EUR,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_EUR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_EUR_boot <- boot_beta$t
    beta_CV_se_validation_raw_EUR <- sd(boot_beta$t)
    
    beta_RV_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_EUR,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_EUR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_EUR_boot <- boot_beta$t
    beta_RV_se_validation_raw_EUR <- sd(boot_beta$t)
    
    auc_validation_raw_EUR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_raw_EUR[!is.na(CV_RV_PRS_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_raw_EUR, statistic = AUC_Boot, R = 10000)
    AUC_raw_EUR_boot <- boot_auc$t
    auc_se_validation_raw_EUR <- sd(boot_auc$t)
    
    beta_CV_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_SAS,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_SAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_SAS_boot <- boot_beta$t
    beta_CV_se_validation_raw_SAS <- sd(boot_beta$t)
    
    beta_RV_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_SAS,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_SAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_SAS_boot <- boot_beta$t
    beta_RV_se_validation_raw_SAS <- sd(boot_beta$t)
    
    auc_validation_raw_SAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_raw_SAS[!is.na(CV_RV_PRS_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_raw_SAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_SAS_boot <- boot_auc$t
    auc_se_validation_raw_SAS <- sd(boot_auc$t)
    
    beta_CV_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_AMR,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_AMR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_AMR_boot <- boot_beta$t
    beta_CV_se_validation_raw_AMR <- sd(boot_beta$t)
    
    beta_RV_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_AMR,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_AMR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_AMR_boot <- boot_beta$t
    beta_RV_se_validation_raw_AMR <- sd(boot_beta$t)
    
    auc_validation_raw_AMR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_raw_AMR[!is.na(CV_RV_PRS_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_raw_AMR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AMR_boot <- boot_auc$t
    auc_se_validation_raw_AMR <- sd(boot_auc$t)
    
    beta_CV_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_AFR,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_AFR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_AFR_boot <- boot_beta$t
    beta_CV_se_validation_raw_AFR <- sd(boot_beta$t)
    
    beta_RV_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_AFR,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_AFR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_AFR_boot <- boot_beta$t
    beta_RV_se_validation_raw_AFR <- sd(boot_beta$t)
    
    auc_validation_raw_AFR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_raw_AFR[!is.na(CV_RV_PRS_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_raw_AFR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AFR_boot <- boot_auc$t
    auc_se_validation_raw_AFR <- sd(boot_auc$t)
    
    beta_CV_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_EAS,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_raw_EAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_raw_EAS_boot <- boot_beta$t
    beta_CV_se_validation_raw_EAS <- sd(boot_beta$t)
    
    beta_RV_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_raw_EAS,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_raw_EAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_raw_EAS_boot <- boot_beta$t
    beta_RV_se_validation_raw_EAS <- sd(boot_beta$t)
    
    auc_validation_raw_EAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_raw_EAS[!is.na(CV_RV_PRS_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_raw_EAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_EAS_boot <- boot_auc$t
    auc_se_validation_raw_EAS <- sd(boot_auc$t)
    
    beta_CV_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_EUR,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_EUR_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_EUR,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_EUR_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    auc_validation_adjusted_EUR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_adjusted_EUR[!is.na(CV_RV_PRS_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EUR_boot <- boot_auc$t
    auc_se_validation_adjusted_EUR <- sd(boot_auc$t)
    
    beta_CV_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_SAS,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_SAS_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_SAS,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_SAS_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    auc_validation_adjusted_SAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_adjusted_SAS[!is.na(CV_RV_PRS_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_SAS_boot <- boot_auc$t
    auc_se_validation_adjusted_SAS <- sd(boot_auc$t)
    
    beta_CV_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_AMR,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_AMR_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_AMR,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_AMR_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    auc_validation_adjusted_AMR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_adjusted_AMR[!is.na(CV_RV_PRS_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AMR_boot <- boot_auc$t
    auc_se_validation_adjusted_AMR <- sd(boot_auc$t)
    
    beta_CV_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_AFR,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_AFR_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_AFR,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_AFR_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    auc_validation_adjusted_AFR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_adjusted_AFR[!is.na(CV_RV_PRS_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AFR_boot <- boot_auc$t
    auc_se_validation_adjusted_AFR <- sd(boot_auc$t)
    
    beta_CV_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_EAS,family = binomial()))[2]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = Beta_CV_Boot, R = 10000)
    beta_CV_adjusted_EAS_boot <- boot_beta$t
    beta_CV_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    beta_RV_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~","CV_PRS + RV_PRS","+",gsub("~","",confounders))),data = CV_RV_PRS_adjusted_EAS,family = binomial()))[3]
    boot_beta <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = Beta_RV_Boot, R = 10000)
    beta_RV_adjusted_EAS_boot <- boot_beta$t
    beta_RV_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    auc_validation_adjusted_EAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = CV_RV_PRS_adjusted_EAS[!is.na(CV_RV_PRS_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EAS_boot <- boot_auc$t
    auc_se_validation_adjusted_EAS <- sd(boot_auc$t)
    
    CV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                 beta_raw = c(beta_CV_validation_raw_EUR,beta_CV_validation_raw_SAS,beta_CV_validation_raw_AMR,beta_CV_validation_raw_AFR,beta_CV_validation_raw_EAS), 
                                 beta_se_raw = c(beta_CV_se_validation_raw_EUR,beta_CV_se_validation_raw_SAS,beta_CV_se_validation_raw_AMR,beta_CV_se_validation_raw_AFR,beta_CV_se_validation_raw_EAS), 
                                 AUC_raw = c(auc_validation_raw_EUR,auc_validation_raw_SAS,auc_validation_raw_AMR,auc_validation_raw_AFR,auc_validation_raw_EAS),
                                 AUC_se_raw = c(auc_se_validation_raw_EUR,auc_se_validation_raw_SAS,auc_se_validation_raw_AMR,auc_se_validation_raw_AFR,auc_se_validation_raw_EAS),
                                 beta_adjusted = c(beta_CV_validation_adjusted_EUR,beta_CV_validation_adjusted_SAS,beta_CV_validation_adjusted_AMR,beta_CV_validation_adjusted_AFR,beta_CV_validation_adjusted_EAS), 
                                 beta_se_adjusted = c(beta_CV_se_validation_adjusted_EUR,beta_CV_se_validation_adjusted_SAS,beta_CV_se_validation_adjusted_AMR,beta_CV_se_validation_adjusted_AFR,beta_CV_se_validation_adjusted_EAS), 
                                 AUC_adjusted = c(auc_validation_adjusted_EUR,auc_validation_adjusted_SAS,auc_validation_adjusted_AMR,auc_validation_adjusted_AFR,auc_validation_adjusted_EAS),
                                 AUC_se_adjusted = c(auc_se_validation_adjusted_EUR,auc_se_validation_adjusted_SAS,auc_se_validation_adjusted_AMR,auc_se_validation_adjusted_AFR,auc_se_validation_adjusted_EAS))
    
    CV_Boot_Results <- data.frame(trait = trait,beta_CV_raw_EUR_boot,AUC_raw_EUR_boot,beta_CV_raw_SAS_boot,AUC_raw_SAS_boot,
                                  beta_CV_raw_AMR_boot,AUC_raw_AMR_boot,beta_CV_raw_AFR_boot,AUC_raw_AFR_boot,
                                  beta_CV_raw_EAS_boot,AUC_raw_EAS_boot,beta_CV_adjusted_EUR_boot,AUC_adjusted_EUR_boot,
                                  beta_CV_adjusted_SAS_boot,AUC_adjusted_SAS_boot,beta_CV_adjusted_AMR_boot,AUC_adjusted_AMR_boot,
                                  beta_CV_adjusted_AFR_boot,AUC_adjusted_AFR_boot,beta_CV_adjusted_EAS_boot,AUC_adjusted_EAS_boot)
    
    RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                 beta_raw = c(beta_RV_validation_raw_EUR,beta_RV_validation_raw_SAS,beta_RV_validation_raw_AMR,beta_RV_validation_raw_AFR,beta_RV_validation_raw_EAS), 
                                 beta_se_raw = c(beta_RV_se_validation_raw_EUR,beta_RV_se_validation_raw_SAS,beta_RV_se_validation_raw_AMR,beta_RV_se_validation_raw_AFR,beta_RV_se_validation_raw_EAS), 
                                 AUC_raw = c(auc_validation_raw_EUR,auc_validation_raw_SAS,auc_validation_raw_AMR,auc_validation_raw_AFR,auc_validation_raw_EAS),
                                 AUC_se_raw = c(auc_se_validation_raw_EUR,auc_se_validation_raw_SAS,auc_se_validation_raw_AMR,auc_se_validation_raw_AFR,auc_se_validation_raw_EAS),
                                 beta_adjusted = c(beta_RV_validation_adjusted_EUR,beta_RV_validation_adjusted_SAS,beta_RV_validation_adjusted_AMR,beta_RV_validation_adjusted_AFR,beta_RV_validation_adjusted_EAS), 
                                 beta_se_adjusted = c(beta_RV_se_validation_adjusted_EUR,beta_RV_se_validation_adjusted_SAS,beta_RV_se_validation_adjusted_AMR,beta_RV_se_validation_adjusted_AFR,beta_RV_se_validation_adjusted_EAS), 
                                 AUC_adjusted = c(auc_validation_adjusted_EUR,auc_validation_adjusted_SAS,auc_validation_adjusted_AMR,auc_validation_adjusted_AFR,auc_validation_adjusted_EAS),
                                 AUC_se_adjusted = c(auc_se_validation_adjusted_EUR,auc_se_validation_adjusted_SAS,auc_se_validation_adjusted_AMR,auc_se_validation_adjusted_AFR,auc_se_validation_adjusted_EAS))
    
    RV_Boot_Results <- data.frame(trait = trait,beta_RV_raw_EUR_boot,AUC_raw_EUR_boot,beta_RV_raw_SAS_boot,AUC_raw_SAS_boot,
                                  beta_RV_raw_AMR_boot,AUC_raw_AMR_boot,beta_RV_raw_AFR_boot,AUC_raw_AFR_boot,
                                  beta_RV_raw_EAS_boot,AUC_raw_EAS_boot,beta_RV_adjusted_EUR_boot,AUC_adjusted_EUR_boot,
                                  beta_RV_adjusted_SAS_boot,AUC_adjusted_SAS_boot,beta_RV_adjusted_AMR_boot,AUC_adjusted_AMR_boot,
                                  beta_RV_adjusted_AFR_boot,AUC_adjusted_AFR_boot,beta_RV_adjusted_EAS_boot,AUC_adjusted_EAS_boot)
  }
  write.csv(CV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"),row.names = FALSE)
  write.csv(CV_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"_Bootstraps.csv"),row.names = FALSE)
  write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"),row.names = FALSE)
  write.csv(RV_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"_Bootstraps.csv"),row.names = FALSE)
  
})[3]
save(time, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/",trait,"_Time.RData"))