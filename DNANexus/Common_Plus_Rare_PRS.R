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

# for array in 1 2 3 4 5 6;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/Common_Plus_Rare_PRS.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/Common_Plus_Rare_PRS.sh -icmd="bash Common_Plus_Rare_PRS.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestPRS/ --priority low --instance-type mem1_ssd1_v2_x4
# done

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "BMI"
}else if(trait == 2){
  trait <- "TC"
}else if(trait == 3){
  trait <- "HDL"
}else if(trait == 4){
  trait <- "LDL"
}else if(trait == 5){
  trait <- "logTG"
}else{
  trait <- "Height"
}


RV_PRS <- read.csv(paste0(trait,"_BestPRS.csv"))
system(paste0("rm ",paste0(trait,"_BestPRS.csv")))
CV_PRS <- read.delim(paste0(trait,"_Best_Validation_All.txt"))
system(paste0("rm ",paste0(trait,"_Best_Validation_All.txt")))

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

load("all_phenotypes.RData")
system("rm all_phenotypes.RData")

CV_RV_PRS_raw_EUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
CV_RV_PRS_raw_NonEUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
CV_RV_PRS_raw_UNK <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
CV_RV_PRS_raw_SAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
CV_RV_PRS_raw_MIX <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
CV_RV_PRS_raw_AFR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
CV_RV_PRS_raw_EAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
CV_RV_PRS_adjusted_NonEUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
CV_RV_PRS_adjusted_UNK <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
CV_RV_PRS_adjusted_MIX <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
CV_RV_PRS_adjusted_EAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

CV_RV_PRS_raw_EUR$Y <- scale(CV_RV_PRS_raw_EUR$Y)
CV_RV_PRS_raw_NonEUR$Y <- scale(CV_RV_PRS_raw_NonEUR$Y)
CV_RV_PRS_raw_UNK$Y <- scale(CV_RV_PRS_raw_UNK$Y)
CV_RV_PRS_raw_SAS$Y <- scale(CV_RV_PRS_raw_SAS$Y)
CV_RV_PRS_raw_MIX$Y <- scale(CV_RV_PRS_raw_MIX$Y)
CV_RV_PRS_raw_AFR$Y <- scale(CV_RV_PRS_raw_AFR$Y)
CV_RV_PRS_raw_EAS$Y <- scale(CV_RV_PRS_raw_EAS$Y)

CV_RV_PRS_adjusted_EUR$Y <- scale(CV_RV_PRS_adjusted_EUR$Y)
CV_RV_PRS_adjusted_NonEUR$Y <- scale(CV_RV_PRS_adjusted_NonEUR$Y)
CV_RV_PRS_adjusted_UNK$Y <- scale(CV_RV_PRS_adjusted_UNK$Y)
CV_RV_PRS_adjusted_SAS$Y <- scale(CV_RV_PRS_adjusted_SAS$Y)
CV_RV_PRS_adjusted_MIX$Y <- scale(CV_RV_PRS_adjusted_MIX$Y)
CV_RV_PRS_adjusted_AFR$Y <- scale(CV_RV_PRS_adjusted_AFR$Y)
CV_RV_PRS_adjusted_EAS$Y <- scale(CV_RV_PRS_adjusted_EAS$Y)

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



best_beta_adjusted_CV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[2]
se_beta_adjusted_CV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))$coefficients[2,2]
best_beta_adjusted_RV_EUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))[3]
se_beta_adjusted_RV_EUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EUR))$coefficients[3,2]

best_beta_adjusted_CV_NonEUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_NonEUR))[2]
se_beta_adjusted_CV_NonEUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_NonEUR))$coefficients[2,2]
best_beta_adjusted_RV_NonEUR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_NonEUR))[3]
se_beta_adjusted_RV_NonEUR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_NonEUR))$coefficients[3,2]

best_beta_adjusted_CV_UNK <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_UNK))[2]
se_beta_adjusted_CV_UNK <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_UNK))$coefficients[2,2]
best_beta_adjusted_RV_UNK <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_UNK))[3]
se_beta_adjusted_RV_UNK <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_UNK))$coefficients[3,2]

best_beta_adjusted_CV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[2]
se_beta_adjusted_CV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))[3]
se_beta_adjusted_RV_SAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_SAS))$coefficients[3,2]

best_beta_adjusted_CV_MIX <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_MIX))[2]
se_beta_adjusted_CV_MIX <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_MIX))$coefficients[2,2]
best_beta_adjusted_RV_MIX <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_MIX))[3]
se_beta_adjusted_RV_MIX <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_MIX))$coefficients[3,2]

best_beta_adjusted_CV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[2]
se_beta_adjusted_CV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))[3]
se_beta_adjusted_RV_AFR <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_AFR))$coefficients[3,2]

best_beta_adjusted_CV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))[2]
se_beta_adjusted_CV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))[3]
se_beta_adjusted_RV_EAS <- summary(lm(Y~prs + RV_PRS,data = CV_RV_PRS_adjusted_EAS))$coefficients[3,2]


CV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_CV_EUR,best_beta_raw_CV_NonEUR,best_beta_raw_CV_UNK,best_beta_raw_CV_SAS,best_beta_raw_CV_MIX,best_beta_raw_CV_AFR,best_beta_raw_CV_EAS), 
                             se_raw = c(se_beta_raw_CV_EUR,se_beta_raw_CV_NonEUR,se_beta_raw_CV_UNK,se_beta_raw_CV_SAS,se_beta_raw_CV_MIX,se_beta_raw_CV_AFR,se_beta_raw_CV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_CV_EUR,best_beta_adjusted_CV_NonEUR,best_beta_adjusted_CV_UNK,best_beta_adjusted_CV_SAS,best_beta_adjusted_CV_MIX,best_beta_adjusted_CV_AFR,best_beta_adjusted_CV_EAS), 
                             se_adjusted = c(se_beta_adjusted_CV_EUR,se_beta_adjusted_CV_NonEUR,se_beta_adjusted_CV_UNK,se_beta_adjusted_CV_SAS,se_beta_adjusted_CV_MIX,se_beta_adjusted_CV_AFR,se_beta_adjusted_CV_EAS),
                             Method = "CV")

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_NonEUR,best_beta_raw_RV_UNK,best_beta_raw_RV_SAS,best_beta_raw_RV_MIX,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_NonEUR,se_beta_raw_RV_UNK,se_beta_raw_RV_SAS,se_beta_raw_RV_MIX,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_NonEUR,best_beta_adjusted_RV_UNK,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_MIX,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                             se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_NonEUR,se_beta_adjusted_RV_UNK,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_MIX,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS),
                             Method = "RV")


write.csv(rbind(CV_PRS_Results,RV_PRS_Results),file = paste0(trait,"Best_Betas.csv"),row.names = FALSE)
