rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/Common_Plus_Rare_PRS_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/Common_Plus_Rare_PRS_Binary.sh -icmd="bash Common_Plus_Rare_PRS_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestPRS/ --priority low --instance-type mem1_ssd1_v2_x4
# done

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


if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <-  paste0("age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}


best_beta_raw_CV_EUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))[2]
se_beta_raw_CV_EUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_EUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))[3]
se_beta_raw_RV_EUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_EUR,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_NonEUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_NonEUR,family = binomial()))[2]
se_beta_raw_CV_NonEUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_NonEUR,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_NonEUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_NonEUR,family = binomial()))[3]
se_beta_raw_RV_NonEUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_NonEUR,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_UNK <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_UNK,family = binomial()))[2]
se_beta_raw_CV_UNK <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_UNK,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_UNK <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_UNK,family = binomial()))[3]
se_beta_raw_RV_UNK <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_UNK,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))[2]
se_beta_raw_CV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))[3]
se_beta_raw_RV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_SAS,family = binomial()))$coefficients[3,2]

best_beta_raw_CV_MIX <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_MIX,family = binomial()))[2]
se_beta_raw_CV_MIX <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_MIX,family = binomial()))$coefficients[2,2]
best_beta_raw_RV_MIX <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_MIX,family = binomial()))[3]
se_beta_raw_RV_MIX <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_raw_MIX,family = binomial()))$coefficients[3,2]

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

best_beta_adjusted_CV_NonEUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_NonEUR,family = binomial()))[2]
se_beta_adjusted_CV_NonEUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_NonEUR,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_NonEUR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_NonEUR,family = binomial()))[3]
se_beta_adjusted_RV_NonEUR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_NonEUR,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_UNK <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_UNK,family = binomial()))[2]
se_beta_adjusted_CV_UNK <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_UNK,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_UNK <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_UNK,family = binomial()))[3]
se_beta_adjusted_RV_UNK <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_UNK,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))[2]
se_beta_adjusted_CV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))[3]
se_beta_adjusted_RV_SAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_SAS,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_MIX <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_MIX,family = binomial()))[2]
se_beta_adjusted_CV_MIX <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_MIX,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_MIX <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_MIX,family = binomial()))[3]
se_beta_adjusted_RV_MIX <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_MIX,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_AFR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))[2]
se_beta_adjusted_CV_AFR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))[3]
se_beta_adjusted_RV_AFR <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_AFR,family = binomial()))$coefficients[3,2]

best_beta_adjusted_CV_EAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))[2]
se_beta_adjusted_CV_EAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))[3]
se_beta_adjusted_RV_EAS <- summary(glm(as.formula(paste0("Y~prs + RV_PRS + ",confounders)),data = CV_RV_PRS_adjusted_EAS,family = binomial()))$coefficients[3,2]


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
