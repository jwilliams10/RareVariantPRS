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

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

RV_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/BestPRS",i,".csv"))
CV_PRS <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/BestPRS",i,".csv")))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Validation_All",i,".txt")))

CV_RV_PRS <- inner_join(RV_PRS,CV_PRS)

CV_RV_PRS_raw <- CV_RV_PRS
CV_RV_PRS_adjusted <- CV_RV_PRS

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
CV_RV_PRS_adjusted <- inner_join(CV_RV_PRS_adjusted,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])

for(j in c("RV_PRS","prs")){
  tmp <- data.frame(y = CV_RV_PRS_adjusted[,j],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
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
    CV_RV_PRS_adjusted[,j] <- 0
  }else{
    CV_RV_PRS_adjusted[,j] <- R/sqrt(y_hat)
  }
}

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

CV_RV_PRS_raw_EUR$Y <- scale(CV_RV_PRS_raw_EUR$Y)
CV_RV_PRS_raw_SAS$Y <- scale(CV_RV_PRS_raw_SAS$Y)
CV_RV_PRS_raw_AMR$Y <- scale(CV_RV_PRS_raw_AMR$Y)
CV_RV_PRS_raw_AFR$Y <- scale(CV_RV_PRS_raw_AFR$Y)
CV_RV_PRS_raw_EAS$Y <- scale(CV_RV_PRS_raw_EAS$Y)

CV_RV_PRS_adjusted_EUR$Y <- scale(CV_RV_PRS_adjusted_EUR$Y)
CV_RV_PRS_adjusted_SAS$Y <- scale(CV_RV_PRS_adjusted_SAS$Y)
CV_RV_PRS_adjusted_AMR$Y <- scale(CV_RV_PRS_adjusted_AMR$Y)
CV_RV_PRS_adjusted_AFR$Y <- scale(CV_RV_PRS_adjusted_AFR$Y)
CV_RV_PRS_adjusted_EAS$Y <- scale(CV_RV_PRS_adjusted_EAS$Y)

CV_RV_PRS_raw_EUR$RV_PRS <- scale(CV_RV_PRS_raw_EUR$RV_PRS)
CV_RV_PRS_raw_SAS$RV_PRS <- scale(CV_RV_PRS_raw_SAS$RV_PRS)
CV_RV_PRS_raw_AMR$RV_PRS <- scale(CV_RV_PRS_raw_AMR$RV_PRS)
CV_RV_PRS_raw_AFR$RV_PRS <- scale(CV_RV_PRS_raw_AFR$RV_PRS)
CV_RV_PRS_raw_EAS$RV_PRS <- scale(CV_RV_PRS_raw_EAS$RV_PRS)

CV_RV_PRS_raw_EUR$prs <- scale(CV_RV_PRS_raw_EUR$prs)
CV_RV_PRS_raw_SAS$prs <- scale(CV_RV_PRS_raw_SAS$prs)
CV_RV_PRS_raw_AMR$prs <- scale(CV_RV_PRS_raw_AMR$prs)
CV_RV_PRS_raw_AFR$prs <- scale(CV_RV_PRS_raw_AFR$prs)
CV_RV_PRS_raw_EAS$prs <- scale(CV_RV_PRS_raw_EAS$prs)


best_beta_raw_CV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))[2]
se_beta_raw_CV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))$coefficients[2,2]
best_beta_raw_RV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))[3]
se_beta_raw_RV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EUR))$coefficients[3,2]

best_beta_raw_CV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))[2]
se_beta_raw_CV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))[3]
se_beta_raw_RV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_SAS))$coefficients[3,2]

best_beta_raw_CV_AMR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AMR))[2]
se_beta_raw_CV_AMR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AMR))$coefficients[2,2]
best_beta_raw_RV_AMR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AMR))[3]
se_beta_raw_RV_AMR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AMR))$coefficients[3,2]

best_beta_raw_CV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))[2]
se_beta_raw_CV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))$coefficients[2,2]
best_beta_raw_RV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))[3]
se_beta_raw_RV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_AFR))$coefficients[3,2]

best_beta_raw_CV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))[2]
se_beta_raw_CV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))$coefficients[2,2]
best_beta_raw_RV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))[3]
se_beta_raw_RV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_raw_EAS))$coefficients[3,2]



best_beta_adjusted_CV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[2]
se_beta_adjusted_CV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))$coefficients[2,2]
best_beta_adjusted_RV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[3]
se_beta_adjusted_RV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))$coefficients[3,2]

best_beta_adjusted_CV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[2]
se_beta_adjusted_CV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[3]
se_beta_adjusted_RV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))$coefficients[3,2]

best_beta_adjusted_CV_AMR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AMR))[2]
se_beta_adjusted_CV_AMR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AMR))$coefficients[2,2]
best_beta_adjusted_RV_AMR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AMR))[3]
se_beta_adjusted_RV_AMR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AMR))$coefficients[3,2]

best_beta_adjusted_CV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[2]
se_beta_adjusted_CV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[3]
se_beta_adjusted_RV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))$coefficients[3,2]

best_beta_adjusted_CV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))[2]
se_beta_adjusted_CV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))[3]
se_beta_adjusted_RV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))$coefficients[3,2]


CV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_CV_EUR,best_beta_raw_CV_SAS,best_beta_raw_CV_AMR,best_beta_raw_CV_AFR,best_beta_raw_CV_EAS), 
                             se_raw = c(se_beta_raw_CV_EUR,se_beta_raw_CV_SAS,se_beta_raw_CV_AMR,se_beta_raw_CV_AFR,se_beta_raw_CV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_CV_EUR,best_beta_adjusted_CV_SAS,best_beta_adjusted_CV_AMR,best_beta_adjusted_CV_AFR,best_beta_adjusted_CV_EAS), 
                             se_adjusted = c(se_beta_adjusted_CV_EUR,se_beta_adjusted_CV_SAS,se_beta_adjusted_CV_AMR,se_beta_adjusted_CV_AFR,se_beta_adjusted_CV_EAS),
                             Method = "CV")

RV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_SAS,best_beta_raw_RV_AMR,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_SAS,se_beta_raw_RV_AMR,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_AMR,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                             se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_AMR,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS),
                             Method = "RV")


write.csv(rbind(CV_PRS_Results,RV_PRS_Results),file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Best_Betas",i,".csv"),row.names = FALSE)
