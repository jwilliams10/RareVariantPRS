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

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

RV_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/BestPRS",i,".csv"))
CV_PRS <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/BestPRS",i,".csv")))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Validation_All",i,".txt")))

CV_RV_PRS <- inner_join(RV_PRS,CV_PRS)

CV_RV_PRS_raw <- CV_RV_PRS

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

CV_RV_PRS_raw_EUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
CV_RV_PRS_raw_NonEUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
CV_RV_PRS_raw_UNK <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
CV_RV_PRS_raw_SAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
CV_RV_PRS_raw_MIX <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
CV_RV_PRS_raw_AFR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
CV_RV_PRS_raw_EAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

CV_RV_PRS_raw_EUR$Y <- scale(CV_RV_PRS_raw_EUR$Y)
CV_RV_PRS_raw_NonEUR$Y <- scale(CV_RV_PRS_raw_NonEUR$Y)
CV_RV_PRS_raw_UNK$Y <- scale(CV_RV_PRS_raw_UNK$Y)
CV_RV_PRS_raw_SAS$Y <- scale(CV_RV_PRS_raw_SAS$Y)
CV_RV_PRS_raw_MIX$Y <- scale(CV_RV_PRS_raw_MIX$Y)
CV_RV_PRS_raw_AFR$Y <- scale(CV_RV_PRS_raw_AFR$Y)
CV_RV_PRS_raw_EAS$Y <- scale(CV_RV_PRS_raw_EAS$Y)

CV_RV_PRS_raw_EUR$RV_PRS <- scale(CV_RV_PRS_raw_EUR$RV_PRS)
CV_RV_PRS_raw_NonEUR$RV_PRS <- scale(CV_RV_PRS_raw_NonEUR$RV_PRS)
CV_RV_PRS_raw_UNK$RV_PRS <- scale(CV_RV_PRS_raw_UNK$RV_PRS)
CV_RV_PRS_raw_SAS$RV_PRS <- scale(CV_RV_PRS_raw_SAS$RV_PRS)
CV_RV_PRS_raw_MIX$RV_PRS <- scale(CV_RV_PRS_raw_MIX$RV_PRS)
CV_RV_PRS_raw_AFR$RV_PRS <- scale(CV_RV_PRS_raw_AFR$RV_PRS)
CV_RV_PRS_raw_EAS$RV_PRS <- scale(CV_RV_PRS_raw_EAS$RV_PRS)

CV_RV_PRS_raw_EUR$prs <- scale(CV_RV_PRS_raw_EUR$prs)
CV_RV_PRS_raw_NonEUR$prs <- scale(CV_RV_PRS_raw_NonEUR$prs)
CV_RV_PRS_raw_UNK$prs <- scale(CV_RV_PRS_raw_UNK$prs)
CV_RV_PRS_raw_SAS$prs <- scale(CV_RV_PRS_raw_SAS$prs)
CV_RV_PRS_raw_MIX$prs <- scale(CV_RV_PRS_raw_MIX$prs)
CV_RV_PRS_raw_AFR$prs <- scale(CV_RV_PRS_raw_AFR$prs)
CV_RV_PRS_raw_EAS$prs <- scale(CV_RV_PRS_raw_EAS$prs)


best_beta_raw_CV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))[2]
se_beta_raw_CV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))$coefficients[2,2]
best_beta_raw_RV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))[3]
se_beta_raw_RV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))$coefficients[3,2]

best_beta_raw_CV_NonEUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_NonEUR))[2]
se_beta_raw_CV_NonEUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_NonEUR))$coefficients[2,2]
best_beta_raw_RV_NonEUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_NonEUR))[3]
se_beta_raw_RV_NonEUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_NonEUR))$coefficients[3,2]

best_beta_raw_CV_UNK <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_UNK))[2]
se_beta_raw_CV_UNK <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_UNK))$coefficients[2,2]
best_beta_raw_RV_UNK <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_UNK))[3]
se_beta_raw_RV_UNK <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_UNK))$coefficients[3,2]

best_beta_raw_CV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))[2]
se_beta_raw_CV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))[3]
se_beta_raw_RV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))$coefficients[3,2]

best_beta_raw_CV_MIX <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_MIX))[2]
se_beta_raw_CV_MIX <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_MIX))$coefficients[2,2]
best_beta_raw_RV_MIX <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_MIX))[3]
se_beta_raw_RV_MIX <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_MIX))$coefficients[3,2]

best_beta_raw_CV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))[2]
se_beta_raw_CV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))$coefficients[2,2]
best_beta_raw_RV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))[3]
se_beta_raw_RV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))$coefficients[3,2]

best_beta_raw_CV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))[2]
se_beta_raw_CV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))$coefficients[2,2]
best_beta_raw_RV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))[3]
se_beta_raw_RV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))$coefficients[3,2]



CV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_CV_EUR,best_beta_raw_CV_NonEUR,best_beta_raw_CV_UNK,best_beta_raw_CV_SAS,best_beta_raw_CV_MIX,best_beta_raw_CV_AFR,best_beta_raw_CV_EAS), 
                             se_raw = c(se_beta_raw_CV_EUR,se_beta_raw_CV_NonEUR,se_beta_raw_CV_UNK,se_beta_raw_CV_SAS,se_beta_raw_CV_MIX,se_beta_raw_CV_AFR,se_beta_raw_CV_EAS),
                             Method = "CV")

RV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_NonEUR,best_beta_raw_RV_UNK,best_beta_raw_RV_SAS,best_beta_raw_RV_MIX,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_NonEUR,se_beta_raw_RV_UNK,se_beta_raw_RV_SAS,se_beta_raw_RV_MIX,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS),
                             Method = "RV")


write.csv(rbind(CV_PRS_Results,RV_PRS_Results),file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Best_Betas",i,".csv"),row.names = FALSE)
