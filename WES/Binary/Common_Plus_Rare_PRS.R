rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "Asthma"
}else if(trait == 2){
  trait <- "CAD"
}else if(trait == 3){
  trait <- "T2D"
}else if(trait == 4){
  trait <- "Breast"
}else{
  trait <- "Prostate"
}

RV_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_BestPRS.csv"))
CV_PRS <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))

CV_RV_PRS <- inner_join(RV_PRS,CV_PRS)

CV_RV_PRS_raw <- CV_RV_PRS
CV_RV_PRS_adjusted <- CV_RV_PRS

for(i in c("RV_PRS","prs")){
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

CV_RV_PRS_raw_EUR$prs <- scale(CV_RV_PRS_raw_EUR$prs)
CV_RV_PRS_raw_SAS$prs <- scale(CV_RV_PRS_raw_SAS$prs)
CV_RV_PRS_raw_AMR$prs <- scale(CV_RV_PRS_raw_AMR$prs)
CV_RV_PRS_raw_AFR$prs <- scale(CV_RV_PRS_raw_AFR$prs)
CV_RV_PRS_raw_EAS$prs <- scale(CV_RV_PRS_raw_EAS$prs)


if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <-  paste0("age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}


best_beta_raw_CV_EUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))[2]
se_beta_raw_CV_EUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_EUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))[3]
se_beta_raw_RV_EUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))[2]
se_beta_raw_CV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))[3]
se_beta_raw_RV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_AMR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AMR,family = binomial()))[2]
se_beta_raw_CV_AMR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AMR,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_AMR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AMR,family = binomial()))[3]
se_beta_raw_RV_AMR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AMR,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_AFR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AFR,family = binomial()))[2]
se_beta_raw_CV_AFR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AFR,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_AFR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AFR,family = binomial()))[3]
se_beta_raw_RV_AFR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_AFR,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_EAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EAS,family = binomial()))[2]
se_beta_raw_CV_EAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EAS,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_EAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EAS,family = binomial()))[3]
se_beta_raw_RV_EAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EAS,family = binomial()))$coefficients[3,2]



best_beta_adjusted_CV_EUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EUR,family = binomial()))[2]
se_beta_adjusted_CV_EUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EUR,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_EUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EUR,family = binomial()))[3]
se_beta_adjusted_RV_EUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EUR,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))[2]
se_beta_adjusted_CV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))[3]
se_beta_adjusted_RV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_AMR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AMR,family = binomial()))[2]
se_beta_adjusted_CV_AMR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AMR,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_AMR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AMR,family = binomial()))[3]
se_beta_adjusted_RV_AMR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AMR,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_AFR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))[2]
se_beta_adjusted_CV_AFR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))[3]
se_beta_adjusted_RV_AFR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_EAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))[2]
se_beta_adjusted_CV_EAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))[3]
se_beta_adjusted_RV_EAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))$coefficients[3,2]


CV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_CV_EUR,best_beta_raw_CV_SAS,best_beta_raw_CV_AMR,best_beta_raw_CV_AFR,best_beta_raw_CV_EAS), 
                             se_raw = c(se_beta_raw_CV_EUR,se_beta_raw_CV_SAS,se_beta_raw_CV_AMR,se_beta_raw_CV_AFR,se_beta_raw_CV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_CV_EUR,best_beta_adjusted_CV_SAS,best_beta_adjusted_CV_AMR,best_beta_adjusted_CV_AFR,best_beta_adjusted_CV_EAS), 
                             se_adjusted = c(se_beta_adjusted_CV_EUR,se_beta_adjusted_CV_SAS,se_beta_adjusted_CV_AMR,se_beta_adjusted_CV_AFR,se_beta_adjusted_CV_EAS),
                             Method = "CV")

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_SAS,best_beta_raw_RV_AMR,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_SAS,se_beta_raw_RV_AMR,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_AMR,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                             se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_AMR,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS),
                             Method = "RV")


write.csv(rbind(CV_PRS_Results,RV_PRS_Results),file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"Best_Betas.csv"),row.names = FALSE)
