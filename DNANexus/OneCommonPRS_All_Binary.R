rm(list = ls())

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/CommonVariant_PRS/OneCommonPRS_All_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/CommonVariant_PRS/OneCommonPRS_All_Binary.sh -icmd="bash OneCommonPRS_All_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/ --priority low --instance-type mem1_ssd1_v2_x4
# done

if(!("caret" %in% rownames(installed.packages()))){
  install.packages("caret",quiet = TRUE)
}

if(!("ranger" %in% rownames(installed.packages()))){
  install.packages("ranger",quiet = TRUE)
}

if(!("RISCA" %in% rownames(installed.packages()))){
  install.packages("RISCA",quiet = TRUE)
}

if(!("SuperLearner" %in% rownames(installed.packages()))){
  install.packages("SuperLearner",quiet = TRUE)
}

if(!("dplyr" %in% rownames(installed.packages()))){
  install.packages("dplyr",quiet = TRUE)
}

if(!("boot" %in% rownames(installed.packages()))){
  install.packages("boot",quiet = TRUE)
}

if(!("stringr" %in% rownames(installed.packages()))){
  install.packages("stringr",quiet = TRUE)
}

library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)
library(stringr)
library(glmnet)

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

#Train
prs_train_CT <- read.csv(paste0(trait,"_prs_all_train.txt"), sep="")
prs_train_LDPred2 <- read.delim(paste0(trait,"_prs_train_ldpred2.sscore"))
prs_train_LASSOSum2 <- read.delim(paste0(trait,"_prs_train_lassosum2.sscore"))

file.remove(paste0(trait,"_prs_all_train.txt"))
file.remove(paste0(trait,"_prs_train_ldpred2.sscore"))
file.remove(paste0(trait,"_prs_train_lassosum2.sscore"))

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],(-1)*prs_train_LDPred2[,-c(1,2,3,4)],(-1)*prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_train <- read.delim("All_Train.txt")
file.remove("All_Train.txt")
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[!is.na(pheno_train[,trait]),-c(1:27)]
pheno_train <- pheno_train[!is.na(pheno_train[,trait]),]

#Tune
prs_tune_CT <- read.csv(paste0(trait,"_prs_all_tune.txt"), sep="")
prs_tune_LDPred2 <- read.delim(paste0(trait,"_prs_tune_ldpred2.sscore"))
prs_tune_LASSOSum2 <- read.delim(paste0(trait,"_prs_tune_lassosum2.sscore"))

file.remove(paste0(trait,"_prs_all_tune.txt"))
file.remove(paste0(trait,"_prs_tune_ldpred2.sscore"))
file.remove(paste0(trait,"_prs_tune_lassosum2.sscore"))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],(-1)*prs_tune_LDPred2[,-c(1,2,3,4)],(-1)*prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_tune <- read.delim("All_Tune.txt")
file.remove("All_Tune.txt")
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[!is.na(pheno_tune[,trait]),-c(1:27)]
pheno_tune <- pheno_tune[!is.na(pheno_tune[,trait]),]

#Validation
prs_validation_CT <- read.csv(paste0(trait,"_prs_all_validation.txt"), sep="")
prs_validation_LDPred2 <- read.delim(paste0(trait,"_prs_validation_ldpred2.sscore"))
prs_validation_LASSOSum2 <- read.delim(paste0(trait,"_prs_validation_lassosum2.sscore"))

file.remove(paste0(trait,"_prs_all_validation.txt"))
file.remove(paste0(trait,"_prs_validation_ldpred2.sscore"))
file.remove(paste0(trait,"_prs_validation_lassosum2.sscore"))

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],(-1)*prs_validation_LDPred2[,-c(1,2,3,4)],(-1)*prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_validation <- read.delim("All_Validation.txt")
file.remove("All_Validation.txt")
pheno_validation <- left_join(pheno_validation,prs_validation_all,by = "IID")

prs_validation_all <- pheno_validation[,-c(1:27)]

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

Results <- Ensemble_Function(x = prs_tune_all,y = pheno_tune[,trait],family = "binary")
Results$Coefficients[is.na(Results$Coefficients)] <- 0
write.csv(Results$Coefficients,file = paste0(trait,"_Final_Coefficients.csv"),row.names = FALSE)
PRS_Train <- as.matrix(pheno_train[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Tune <- as.matrix(pheno_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Validation <- as.matrix(pheno_validation[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)

prs_best_train <- data.frame(IID = pheno_train$IID,prs = PRS_Train)
prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = PRS_Tune)
prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = PRS_Validation)

write.table(prs_best_train,file=paste0(trait,"_Best_Train_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0(trait,"_Best_Tune_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0(trait,"_Best_Validation_All.txt"),sep = "\t",row.names = FALSE)


if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <- paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}

load("all_phenotypes.RData")

pheno_validation_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),trait)],prs_best_validation)
pheno_validation_adjusted <- pheno_validation_raw
tmp <- data.frame(y = pheno_validation_adjusted[,"prs"],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
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
  pheno_validation_adjusted[,"prs"] <- 0
}else{
  pheno_validation_adjusted[,"prs"] <- R/sqrt(y_hat)
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

pheno_validation_raw_EUR[,"prs"] <- scale(pheno_validation_raw_EUR[,"prs"])
pheno_validation_raw_SAS[,"prs"] <- scale(pheno_validation_raw_SAS[,"prs"])
pheno_validation_raw_AMR[,"prs"] <- scale(pheno_validation_raw_AMR[,"prs"])
pheno_validation_raw_AFR[,"prs"] <- scale(pheno_validation_raw_AFR[,"prs"])
pheno_validation_raw_EAS[,"prs"] <- scale(pheno_validation_raw_EAS[,"prs"])


Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
  return(c(result))
}

AUC_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  return(c(result))
}

beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

auc_validation_raw_EUR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_EUR[!is.na(pheno_validation_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_EUR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_EUR <- sd(boot_auc$t)
auc_lower_validation_raw_EUR <- auc_ci$basic[4]
auc_upper_validation_raw_EUR <- auc_ci$basic[5]

beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

auc_validation_raw_SAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_SAS[!is.na(pheno_validation_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_SAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_SAS <- sd(boot_auc$t)
auc_lower_validation_raw_SAS <- auc_ci$basic[4]
auc_upper_validation_raw_SAS <- auc_ci$basic[5]

beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

auc_validation_raw_AMR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_AMR[!is.na(pheno_validation_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_AMR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_AMR <- sd(boot_auc$t)
auc_lower_validation_raw_AMR <- auc_ci$basic[4]
auc_upper_validation_raw_AMR <- auc_ci$basic[5]

beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

auc_validation_raw_AFR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_AFR[!is.na(pheno_validation_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_AFR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_AFR <- sd(boot_auc$t)
auc_lower_validation_raw_AFR <- auc_ci$basic[4]
auc_upper_validation_raw_AFR <- auc_ci$basic[5]

beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EAS <- sd(boot_beta$t)
beta_lower_validation_raw_EAS <- beta_ci$basic[4]
beta_upper_validation_raw_EAS <- beta_ci$basic[5]

auc_validation_raw_EAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_EAS[!is.na(pheno_validation_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_EAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_EAS <- sd(boot_auc$t)
auc_lower_validation_raw_EAS <- auc_ci$basic[4]
auc_upper_validation_raw_EAS <- auc_ci$basic[5]

beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

auc_validation_adjusted_EUR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_EUR[!is.na(pheno_validation_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_EUR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_EUR <- sd(boot_auc$t)
auc_lower_validation_adjusted_EUR <- auc_ci$basic[4]
auc_upper_validation_adjusted_EUR <- auc_ci$basic[5]

beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

auc_validation_adjusted_SAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_SAS[!is.na(pheno_validation_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_SAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_SAS <- sd(boot_auc$t)
auc_lower_validation_adjusted_SAS <- auc_ci$basic[4]
auc_upper_validation_adjusted_SAS <- auc_ci$basic[5]

beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

auc_validation_adjusted_AMR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_AMR[!is.na(pheno_validation_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_AMR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_AMR <- sd(boot_auc$t)
auc_lower_validation_adjusted_AMR <- auc_ci$basic[4]
auc_upper_validation_adjusted_AMR <- auc_ci$basic[5]

beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

auc_validation_adjusted_AFR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_AFR[!is.na(pheno_validation_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_AFR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_AFR <- sd(boot_auc$t)
auc_lower_validation_adjusted_AFR <- auc_ci$basic[4]
auc_upper_validation_adjusted_AFR <- auc_ci$basic[5]

beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]

auc_validation_adjusted_EAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_EAS[!is.na(pheno_validation_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_EAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_EAS <- sd(boot_auc$t)
auc_lower_validation_adjusted_EAS <- auc_ci$basic[4]
auc_upper_validation_adjusted_EAS <- auc_ci$basic[5]

SL_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                         beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                         beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR,beta_lower_validation_raw_EAS), 
                         beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR,beta_upper_validation_raw_EAS), 
                         AUC_raw = c(auc_validation_raw_EUR,auc_validation_raw_SAS,auc_validation_raw_AMR,auc_validation_raw_AFR,auc_validation_raw_EAS),
                         AUC_se_raw = c(auc_se_validation_raw_EUR,auc_se_validation_raw_SAS,auc_se_validation_raw_AMR,auc_se_validation_raw_AFR,auc_se_validation_raw_EAS),
                         AUC_lower_raw = c(auc_lower_validation_raw_EUR,auc_lower_validation_raw_SAS,auc_lower_validation_raw_AMR,auc_lower_validation_raw_AFR,auc_lower_validation_raw_EAS),
                         AUC_upper_raw = c(auc_upper_validation_raw_EUR,auc_upper_validation_raw_SAS,auc_upper_validation_raw_AMR,auc_upper_validation_raw_AFR,auc_upper_validation_raw_EAS),
                         beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                         beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                         beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR,beta_lower_validation_adjusted_EAS), 
                         beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR,beta_upper_validation_adjusted_EAS), 
                         AUC_adjusted = c(auc_validation_adjusted_EUR,auc_validation_adjusted_SAS,auc_validation_adjusted_AMR,auc_validation_adjusted_AFR,auc_validation_adjusted_EAS),
                         AUC_se_adjusted = c(auc_se_validation_adjusted_EUR,auc_se_validation_adjusted_SAS,auc_se_validation_adjusted_AMR,auc_se_validation_adjusted_AFR,auc_se_validation_adjusted_EAS),
                         AUC_lower_adjusted = c(auc_lower_validation_adjusted_EUR,auc_lower_validation_adjusted_SAS,auc_lower_validation_adjusted_AMR,auc_lower_validation_adjusted_AFR,auc_lower_validation_adjusted_EAS),
                         AUC_upper_adjusted = c(auc_upper_validation_adjusted_EUR,auc_upper_validation_adjusted_SAS,auc_upper_validation_adjusted_AMR,auc_upper_validation_adjusted_AFR,auc_upper_validation_adjusted_EAS))

write.csv(SL_Results,file = paste0(trait,"Best_Betas.csv"),row.names = FALSE)