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

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")
pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")
CV_PRS_Tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"))
colnames(CV_PRS_Tune) <- c("IID","CV_PRS")
pheno_tune <- inner_join(pheno_tune,CV_PRS_Tune)
RV_PRS_Tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/PRS_Tune_",i,".csv"))
colnames(RV_PRS_Tune) <- c("IID","RV_PRS")
pheno_tune <- inner_join(pheno_tune,RV_PRS_Tune)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")
pheno_validation <- Y_validation[[i]]
colnames(pheno_validation) <- c("IID","Y")
CV_PRS_Validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))
colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/PRS_Validation_",i,".csv"))
colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Tune_All",i,".txt")))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/PRS_Tune_",i,".csv")))

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Validation_All",i,".txt")))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/PRS_Validation_",i,".csv")))


model.null <- lm(Y~1,data=pheno_tune)
pheno_tune$y_tune <- model.null$residual

model.null <- lm(Y~1,data=pheno_validation)
pheno_validation$y_validation <- model.null$residual

RICE_Model <- lm(y_tune ~ CV_PRS + RV_PRS,data = pheno_tune)
pheno_validation$PRS <- predict(RICE_Model,pheno_validation)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

CV_RV_PRS_raw <- inner_join(pheno_validation,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])
CV_RV_PRS_adjusted <- CV_RV_PRS_raw

for(column_name in c("CV_PRS","RV_PRS","PRS")){
  tmp <- data.frame(y = CV_RV_PRS_adjusted[,column_name],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
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
    CV_RV_PRS_adjusted[,column_name] <- 0
  }else{
    CV_RV_PRS_adjusted[,column_name] <- R/sqrt(y_hat)
  }
}

CV_RV_PRS_raw_EUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
CV_RV_PRS_raw_SAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
CV_RV_PRS_raw_AMR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
CV_RV_PRS_raw_AFR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]

CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]

CV_RV_PRS_raw_EUR$y_validation <- scale(CV_RV_PRS_raw_EUR$y_validation)
CV_RV_PRS_raw_SAS$y_validation <- scale(CV_RV_PRS_raw_SAS$y_validation)
CV_RV_PRS_raw_AMR$y_validation <- scale(CV_RV_PRS_raw_AMR$y_validation)
CV_RV_PRS_raw_AFR$y_validation <- scale(CV_RV_PRS_raw_AFR$y_validation)

CV_RV_PRS_adjusted_EUR$y_validation <- scale(CV_RV_PRS_adjusted_EUR$y_validation)
CV_RV_PRS_adjusted_SAS$y_validation <- scale(CV_RV_PRS_adjusted_SAS$y_validation)
CV_RV_PRS_adjusted_AMR$y_validation <- scale(CV_RV_PRS_adjusted_AMR$y_validation)
CV_RV_PRS_adjusted_AFR$y_validation <- scale(CV_RV_PRS_adjusted_AFR$y_validation)

CV_RV_PRS_raw_EUR$RV_PRS <- scale(CV_RV_PRS_raw_EUR$RV_PRS)
CV_RV_PRS_raw_SAS$RV_PRS <- scale(CV_RV_PRS_raw_SAS$RV_PRS)
CV_RV_PRS_raw_AMR$RV_PRS <- scale(CV_RV_PRS_raw_AMR$RV_PRS)
CV_RV_PRS_raw_AFR$RV_PRS <- scale(CV_RV_PRS_raw_AFR$RV_PRS)

CV_RV_PRS_raw_EUR$CV_PRS <- scale(CV_RV_PRS_raw_EUR$CV_PRS)
CV_RV_PRS_raw_SAS$CV_PRS <- scale(CV_RV_PRS_raw_SAS$CV_PRS)
CV_RV_PRS_raw_AMR$CV_PRS <- scale(CV_RV_PRS_raw_AMR$CV_PRS)
CV_RV_PRS_raw_AFR$CV_PRS <- scale(CV_RV_PRS_raw_AFR$CV_PRS)

CV_RV_PRS_raw_EUR$PRS <- scale(CV_RV_PRS_raw_EUR$PRS)
CV_RV_PRS_raw_SAS$PRS <- scale(CV_RV_PRS_raw_SAS$PRS)
CV_RV_PRS_raw_AMR$PRS <- scale(CV_RV_PRS_raw_AMR$PRS)
CV_RV_PRS_raw_AFR$PRS <- scale(CV_RV_PRS_raw_AFR$PRS)


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
boot_beta <- boot(data = CV_RV_PRS_raw_EUR, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_raw_EUR <- sd(boot_beta$t)
beta_CV_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_CV_upper_validation_raw_EUR <- beta_ci$basic[5]

beta_RV_validation_raw_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_EUR))[3]
boot_beta <- boot(data = CV_RV_PRS_raw_EUR, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_raw_EUR <- sd(boot_beta$t)
beta_RV_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_RV_upper_validation_raw_EUR <- beta_ci$basic[5]

R2_validation_raw_EUR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_EUR))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_raw_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EUR <- sd(boot_R2$t)
R2_lower_validation_raw_EUR <- R2_ci$basic[4]
R2_upper_validation_raw_EUR <- R2_ci$basic[5]

beta_CV_validation_raw_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_SAS))[2]
boot_beta <- boot(data = CV_RV_PRS_raw_SAS, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_raw_SAS <- sd(boot_beta$t)
beta_CV_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_CV_upper_validation_raw_SAS <- beta_ci$basic[5]

beta_RV_validation_raw_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_SAS))[3]
boot_beta <- boot(data = CV_RV_PRS_raw_SAS, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_raw_SAS <- sd(boot_beta$t)
beta_RV_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_RV_upper_validation_raw_SAS <- beta_ci$basic[5]

R2_validation_raw_SAS <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_SAS))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_raw_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_SAS <- sd(boot_R2$t)
R2_lower_validation_raw_SAS <- R2_ci$basic[4]
R2_upper_validation_raw_SAS <- R2_ci$basic[5]

beta_CV_validation_raw_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AMR))[2]
boot_beta <- boot(data = CV_RV_PRS_raw_AMR, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_raw_AMR <- sd(boot_beta$t)
beta_CV_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_CV_upper_validation_raw_AMR <- beta_ci$basic[5]

beta_RV_validation_raw_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AMR))[3]
boot_beta <- boot(data = CV_RV_PRS_raw_AMR, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_raw_AMR <- sd(boot_beta$t)
beta_RV_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_RV_upper_validation_raw_AMR <- beta_ci$basic[5]

R2_validation_raw_AMR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_AMR))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_raw_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AMR <- sd(boot_R2$t)
R2_lower_validation_raw_AMR <- R2_ci$basic[4]
R2_upper_validation_raw_AMR <- R2_ci$basic[5]

beta_CV_validation_raw_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AFR))[2]
boot_beta <- boot(data = CV_RV_PRS_raw_AFR, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_raw_AFR <- sd(boot_beta$t)
beta_CV_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_CV_upper_validation_raw_AFR <- beta_ci$basic[5]

beta_RV_validation_raw_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_raw_AFR))[3]
boot_beta <- boot(data = CV_RV_PRS_raw_AFR, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_raw_AFR <- sd(boot_beta$t)
beta_RV_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_RV_upper_validation_raw_AFR <- beta_ci$basic[5]

R2_validation_raw_AFR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_raw_AFR))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_raw_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AFR <- sd(boot_R2$t)
R2_lower_validation_raw_AFR <- R2_ci$basic[4]
R2_upper_validation_raw_AFR <- R2_ci$basic[5]

beta_CV_validation_adjusted_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[2]
boot_beta <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_CV_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_CV_upper_validation_adjusted_EUR <- beta_ci$basic[5]

beta_RV_validation_adjusted_EUR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[3]
boot_beta <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_RV_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_RV_upper_validation_adjusted_EUR <- beta_ci$basic[5]

R2_validation_adjusted_EUR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_EUR))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
R2_lower_validation_adjusted_EUR <- R2_ci$basic[4]
R2_upper_validation_adjusted_EUR <- R2_ci$basic[5]

beta_CV_validation_adjusted_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[2]
boot_beta <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_CV_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_CV_upper_validation_adjusted_SAS <- beta_ci$basic[5]

beta_RV_validation_adjusted_SAS <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[3]
boot_beta <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_RV_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_RV_upper_validation_adjusted_SAS <- beta_ci$basic[5]

R2_validation_adjusted_SAS <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_SAS))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_SAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_SAS <- R2_ci$basic[5]

beta_CV_validation_adjusted_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AMR))[2]
boot_beta <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_CV_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_CV_upper_validation_adjusted_AMR <- beta_ci$basic[5]

beta_RV_validation_adjusted_AMR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AMR))[3]
boot_beta <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_RV_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_RV_upper_validation_adjusted_AMR <- beta_ci$basic[5]

R2_validation_adjusted_AMR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_AMR))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AMR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AMR <- R2_ci$basic[5]

beta_CV_validation_adjusted_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[2]
boot_beta <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = Beta_CV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_CV_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_CV_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_CV_upper_validation_adjusted_AFR <- beta_ci$basic[5]

beta_RV_validation_adjusted_AFR <- coef(lm(y_validation~CV_PRS + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[3]
boot_beta <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = Beta_RV_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_RV_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_RV_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_RV_upper_validation_adjusted_AFR <- beta_ci$basic[5]

R2_validation_adjusted_AFR <- summary(lm(y_validation~PRS,data = CV_RV_PRS_adjusted_AFR))$r.squared
boot_R2 <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AFR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AFR <- R2_ci$basic[5]

CV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR"), 
                             beta_raw = c(beta_CV_validation_raw_EUR,beta_CV_validation_raw_SAS,beta_CV_validation_raw_AMR,beta_CV_validation_raw_AFR), 
                             beta_se_raw = c(beta_CV_se_validation_raw_EUR,beta_CV_se_validation_raw_SAS,beta_CV_se_validation_raw_AMR,beta_CV_se_validation_raw_AFR), 
                             beta_lower_raw = c(beta_CV_lower_validation_raw_EUR,beta_CV_lower_validation_raw_SAS,beta_CV_lower_validation_raw_AMR,beta_CV_lower_validation_raw_AFR), 
                             beta_upper_raw = c(beta_CV_upper_validation_raw_EUR,beta_CV_upper_validation_raw_SAS,beta_CV_upper_validation_raw_AMR,beta_CV_upper_validation_raw_AFR), 
                             R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR),
                             R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR),
                             R2_lower_raw = c(R2_lower_validation_raw_EUR,R2_lower_validation_raw_SAS,R2_lower_validation_raw_AMR,R2_lower_validation_raw_AFR),
                             R2_upper_raw = c(R2_upper_validation_raw_EUR,R2_upper_validation_raw_SAS,R2_upper_validation_raw_AMR,R2_upper_validation_raw_AFR),
                             beta_adjusted = c(beta_CV_validation_adjusted_EUR,beta_CV_validation_adjusted_SAS,beta_CV_validation_adjusted_AMR,beta_CV_validation_adjusted_AFR), 
                             beta_se_adjusted = c(beta_CV_se_validation_adjusted_EUR,beta_CV_se_validation_adjusted_SAS,beta_CV_se_validation_adjusted_AMR,beta_CV_se_validation_adjusted_AFR), 
                             beta_lower_adjusted = c(beta_CV_lower_validation_adjusted_EUR,beta_CV_lower_validation_adjusted_SAS,beta_CV_lower_validation_adjusted_AMR,beta_CV_lower_validation_adjusted_AFR), 
                             beta_upper_adjusted = c(beta_CV_upper_validation_adjusted_EUR,beta_CV_upper_validation_adjusted_SAS,beta_CV_upper_validation_adjusted_AMR,beta_CV_upper_validation_adjusted_AFR), 
                             R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR),
                             R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR),
                             R2_lower_adjusted = c(R2_lower_validation_adjusted_EUR,R2_lower_validation_adjusted_SAS,R2_lower_validation_adjusted_AMR,R2_lower_validation_adjusted_AFR),
                             R2_upper_adjusted = c(R2_upper_validation_adjusted_EUR,R2_upper_validation_adjusted_SAS,R2_upper_validation_adjusted_AMR,R2_upper_validation_adjusted_AFR))

RV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR"), 
                             beta_raw = c(beta_RV_validation_raw_EUR,beta_RV_validation_raw_SAS,beta_RV_validation_raw_AMR,beta_RV_validation_raw_AFR), 
                             beta_se_raw = c(beta_RV_se_validation_raw_EUR,beta_RV_se_validation_raw_SAS,beta_RV_se_validation_raw_AMR,beta_RV_se_validation_raw_AFR), 
                             beta_lower_raw = c(beta_RV_lower_validation_raw_EUR,beta_RV_lower_validation_raw_SAS,beta_RV_lower_validation_raw_AMR,beta_RV_lower_validation_raw_AFR), 
                             beta_upper_raw = c(beta_RV_upper_validation_raw_EUR,beta_RV_upper_validation_raw_SAS,beta_RV_upper_validation_raw_AMR,beta_RV_upper_validation_raw_AFR), 
                             R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR),
                             R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR),
                             R2_lower_raw = c(R2_lower_validation_raw_EUR,R2_lower_validation_raw_SAS,R2_lower_validation_raw_AMR,R2_lower_validation_raw_AFR),
                             R2_upper_raw = c(R2_upper_validation_raw_EUR,R2_upper_validation_raw_SAS,R2_upper_validation_raw_AMR,R2_upper_validation_raw_AFR),
                             beta_adjusted = c(beta_RV_validation_adjusted_EUR,beta_RV_validation_adjusted_SAS,beta_RV_validation_adjusted_AMR,beta_RV_validation_adjusted_AFR), 
                             beta_se_adjusted = c(beta_RV_se_validation_adjusted_EUR,beta_RV_se_validation_adjusted_SAS,beta_RV_se_validation_adjusted_AMR,beta_RV_se_validation_adjusted_AFR), 
                             beta_lower_adjusted = c(beta_RV_lower_validation_adjusted_EUR,beta_RV_lower_validation_adjusted_SAS,beta_RV_lower_validation_adjusted_AMR,beta_RV_lower_validation_adjusted_AFR), 
                             beta_upper_adjusted = c(beta_RV_upper_validation_adjusted_EUR,beta_RV_upper_validation_adjusted_SAS,beta_RV_upper_validation_adjusted_AMR,beta_RV_upper_validation_adjusted_AFR), 
                             R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR),
                             R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR),
                             R2_lower_adjusted = c(R2_lower_validation_adjusted_EUR,R2_lower_validation_adjusted_SAS,R2_lower_validation_adjusted_AMR,R2_lower_validation_adjusted_AFR),
                             R2_upper_adjusted = c(R2_upper_validation_adjusted_EUR,R2_upper_validation_adjusted_SAS,R2_upper_validation_adjusted_AMR,R2_upper_validation_adjusted_AFR))

write.csv(CV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"),row.names = FALSE)
write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"),row.names = FALSE)
