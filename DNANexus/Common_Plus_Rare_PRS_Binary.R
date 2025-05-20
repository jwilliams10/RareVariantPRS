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

pheno_tune <- read.delim("All_Tune.txt")
CV_PRS_Tune <- read.delim(paste0(trait,"_Best_Tune_All.txt"))
system(paste0("rm ",paste0(trait,"_Best_Tune_All.txt")))
colnames(CV_PRS_Tune) <- c("IID","CV_PRS")
pheno_tune <- inner_join(pheno_tune,CV_PRS_Tune)
RV_PRS_Tune <- read.csv(paste0(trait,"Tune_BestPRS.csv"))
system(paste0("rm ",paste0(trait,"Tune_BestPRS.csv")))
colnames(RV_PRS_Tune) <- c("IID","RV_PRS")
pheno_tune <- inner_join(pheno_tune,RV_PRS_Tune)

pheno_validation <- read.delim("All_Validation.txt")
CV_PRS_Validation <- read.delim(paste0(trait,"_Best_Validation_All.txt"))
system(paste0("rm ",paste0(trait,"_Best_Validation_All.txt")))
colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
RV_PRS_Validation <- read.csv(paste0(trait,"Validation_BestPRS.csv"))
system(paste0("rm ",paste0(trait,"Validation_BestPRS.csv")))
colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
CT_PRS_Validation <- read.delim(paste0(trait,"_prs_validation_best.txt"),sep = " ",header = TRUE)
system(paste0("rm ",paste0(trait,"_prs_validation_best.txt")))
colnames(CT_PRS_Validation) <- c("IID","FID","CT_PRS")
CT_PRS_Validation <- CT_PRS_Validation[,c("IID","CT_PRS")]
pheno_validation <- inner_join(pheno_validation,CT_PRS_Validation)
LDpred2_PRS_Validation <- read.delim(paste0(trait,"_ldpred2_validation_prs_best.txt"),sep = "\t",header = TRUE)
system(paste0("rm ",paste0(trait,"_ldpred2_validation_prs_best.txt")))
colnames(LDpred2_PRS_Validation) <- c("IID","LDpred2_PRS")
LDpred2_PRS_Validation$LDpred2_PRS <- (-1)*LDpred2_PRS_Validation$LDpred2_PRS
pheno_validation <- inner_join(pheno_validation,LDpred2_PRS_Validation)
Lassosum2_PRS_Validation <- read.delim(paste0(trait,"_lassosum2_validation_prs_best.txt"),sep = "\t",header = TRUE)
system(paste0("rm ",paste0(trait,"_lassosum2_validation_prs_best.txt")))
colnames(Lassosum2_PRS_Validation) <- c("IID","Lassosum2_PRS")
Lassosum2_PRS_Validation$Lassosum2_PRS <- (-1)*Lassosum2_PRS_Validation$Lassosum2_PRS
pheno_validation <- inner_join(pheno_validation,Lassosum2_PRS_Validation)

RICE_Model <- glm(as.formula(paste0(trait,"~ CV_PRS + RV_PRS")),data = pheno_tune,family = binomial())
pheno_validation$PRS <- predict(RICE_Model,pheno_validation,type = "link")

CV_RV_PRS_raw <- pheno_validation
CV_RV_PRS_adjusted <- pheno_validation

for(i in c("CV_PRS","RV_PRS","PRS","CT_PRS","LDpred2_PRS","Lassosum2_PRS")){
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

CV_RV_PRS_raw_EUR$CT_PRS <- scale(CV_RV_PRS_raw_EUR$CT_PRS)
CV_RV_PRS_raw_SAS$CT_PRS <- scale(CV_RV_PRS_raw_SAS$CT_PRS)
CV_RV_PRS_raw_AMR$CT_PRS <- scale(CV_RV_PRS_raw_AMR$CT_PRS)
CV_RV_PRS_raw_AFR$CT_PRS <- scale(CV_RV_PRS_raw_AFR$CT_PRS)
CV_RV_PRS_raw_EAS$CT_PRS <- scale(CV_RV_PRS_raw_EAS$CT_PRS)

CV_RV_PRS_raw_EUR$LDpred2_PRS <- scale(CV_RV_PRS_raw_EUR$LDpred2_PRS)
CV_RV_PRS_raw_SAS$LDpred2_PRS <- scale(CV_RV_PRS_raw_SAS$LDpred2_PRS)
CV_RV_PRS_raw_AMR$LDpred2_PRS <- scale(CV_RV_PRS_raw_AMR$LDpred2_PRS)
CV_RV_PRS_raw_AFR$LDpred2_PRS <- scale(CV_RV_PRS_raw_AFR$LDpred2_PRS)
CV_RV_PRS_raw_EAS$LDpred2_PRS <- scale(CV_RV_PRS_raw_EAS$LDpred2_PRS)

CV_RV_PRS_raw_EUR$Lassosum2_PRS <- scale(CV_RV_PRS_raw_EUR$Lassosum2_PRS)
CV_RV_PRS_raw_SAS$Lassosum2_PRS <- scale(CV_RV_PRS_raw_SAS$Lassosum2_PRS)
CV_RV_PRS_raw_AMR$Lassosum2_PRS <- scale(CV_RV_PRS_raw_AMR$Lassosum2_PRS)
CV_RV_PRS_raw_AFR$Lassosum2_PRS <- scale(CV_RV_PRS_raw_AFR$Lassosum2_PRS)
CV_RV_PRS_raw_EAS$Lassosum2_PRS <- scale(CV_RV_PRS_raw_EAS$Lassosum2_PRS)


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

AUC_Comparison_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  RICE_AUC <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  CT_AUC <- roc.binary(status = trait,variable = "CT_PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  LDpred2_AUC <- roc.binary(status = trait,variable = "LDpred2_PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  Lassosum2_AUC <- roc.binary(status = trait,variable = "Lassosum2_PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  return(c(RICE_AUC - CT_AUC,RICE_AUC - LDpred2_AUC,RICE_AUC - Lassosum2_AUC))
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

boot_AUC <- boot(data = CV_RV_PRS_raw_EUR, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_raw_EUR_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_raw_EUR_boot) <- c("AUC_raw_EUR_RICE_vs_CT","AUC_raw_EUR_RICE_vs_LDpred2","AUC_raw_EUR_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_raw_SAS, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_raw_SAS_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_raw_SAS_boot) <- c("AUC_raw_SAS_RICE_vs_CT","AUC_raw_SAS_RICE_vs_LDpred2","AUC_raw_SAS_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_raw_AMR, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_raw_AMR_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_raw_AMR_boot) <- c("AUC_raw_AMR_RICE_vs_CT","AUC_raw_AMR_RICE_vs_LDpred2","AUC_raw_AMR_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_raw_AFR, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_raw_AFR_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_raw_AFR_boot) <- c("AUC_raw_AFR_RICE_vs_CT","AUC_raw_AFR_RICE_vs_LDpred2","AUC_raw_AFR_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_raw_EAS, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_raw_EAS_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_raw_EAS_boot) <- c("AUC_raw_EAS_RICE_vs_CT","AUC_raw_EAS_RICE_vs_LDpred2","AUC_raw_EAS_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_adjusted_EUR, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_adjusted_EUR_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_adjusted_EUR_boot) <- c("AUC_adjusted_EUR_RICE_vs_CT","AUC_adjusted_EUR_RICE_vs_LDpred2","AUC_adjusted_EUR_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_adjusted_SAS, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_adjusted_SAS_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_adjusted_SAS_boot) <- c("AUC_adjusted_SAS_RICE_vs_CT","AUC_adjusted_SAS_RICE_vs_LDpred2","AUC_adjusted_SAS_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_adjusted_AMR, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_adjusted_AMR_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_adjusted_AMR_boot) <- c("AUC_adjusted_AMR_RICE_vs_CT","AUC_adjusted_AMR_RICE_vs_LDpred2","AUC_adjusted_AMR_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_adjusted_AFR, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_adjusted_AFR_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_adjusted_AFR_boot) <- c("AUC_adjusted_AFR_RICE_vs_CT","AUC_adjusted_AFR_RICE_vs_LDpred2","AUC_adjusted_AFR_RICE_vs_Lassosum2")

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

boot_AUC <- boot(data = CV_RV_PRS_adjusted_EAS, statistic = AUC_Comparison_Boot, R = 10000)
AUC_comparison_adjusted_EAS_boot <- as.data.frame(boot_AUC$t)
colnames(AUC_comparison_adjusted_EAS_boot) <- c("AUC_adjusted_EAS_RICE_vs_CT","AUC_adjusted_EAS_RICE_vs_LDpred2","AUC_adjusted_EAS_RICE_vs_Lassosum2")

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

Comparison_Boot_Results <- data.frame(trait = trait,AUC_comparison_raw_EUR_boot,AUC_comparison_raw_SAS_boot,AUC_comparison_raw_AMR_boot,AUC_comparison_raw_AFR_boot,AUC_comparison_raw_EAS_boot,
                                      AUC_comparison_adjusted_EUR_boot,AUC_comparison_adjusted_SAS_boot,AUC_comparison_adjusted_AMR_boot,AUC_comparison_adjusted_AFR_boot,AUC_comparison_adjusted_EAS_boot)

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

write.csv(Comparison_Boot_Results,file = paste0(trait,"_Comparison_Bootstraps.csv"),row.names = FALSE)  
write.csv(CV_PRS_Results,file = paste0("CV_",trait,"Best_Betas.csv"),row.names = FALSE)
write.csv(CV_Boot_Results,file = paste0("CV_",trait,"_Bootstraps.csv"),row.names = FALSE)
write.csv(RV_PRS_Results,file = paste0("RV_",trait,"Best_Betas.csv"),row.names = FALSE)
write.csv(RV_Boot_Results,file = paste0("RV_",trait,"_Bootstraps.csv"),row.names = FALSE)
