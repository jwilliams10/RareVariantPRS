rm(list = ls())

# for array in 1 2 3 4 5 6;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/CommonVariant_PRS/OneCommonPRS_All.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/CommonVariant_PRS/OneCommonPRS_All.sh -icmd="bash OneCommonPRS_All.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/ --priority low --instance-type mem1_ssd1_v2_x4
# done

if(!("caret" %in% rownames(installed.packages()))){
  install.packages("caret",quiet = TRUE)
}

if(!("ranger" %in% rownames(installed.packages()))){
  install.packages("ranger",quiet = TRUE)
}

if(!("bigsparser" %in% rownames(installed.packages()))){
  install.packages("bigsparser",quiet = TRUE)
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
library(stringr)
library(glmnet)

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

#Train
prs_train_CT <- read.csv(paste0(trait,"_prs_all_train.txt"), sep="")
prs_train_LDPred2 <- read.delim(paste0(trait,"_prs_train_ldpred2.sscore"))
prs_train_LASSOSum2 <- read.delim(paste0(trait,"_prs_train_lassosum2.sscore"))

file.remove(paste0(trait,"_prs_all_train.txt"))
file.remove(paste0(trait,"_prs_train_ldpred2.sscore"))
file.remove(paste0(trait,"_prs_train_lassosum2.sscore"))

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],prs_train_LDPred2[,-c(1,2,3,4)],prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_train <- read.delim("All_Train.txt")
file.remove("All_Train.txt")
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[,-c(1:27)]

#Tune
prs_tune_CT <- read.csv(paste0(trait,"_prs_all_tune.txt"), sep="")
prs_tune_LDPred2 <- read.delim(paste0(trait,"_prs_tune_ldpred2.sscore"))
prs_tune_LASSOSum2 <- read.delim(paste0(trait,"_prs_tune_lassosum2.sscore"))

file.remove(paste0(trait,"_prs_all_tune.txt"))
file.remove(paste0(trait,"_prs_tune_ldpred2.sscore"))
file.remove(paste0(trait,"_prs_tune_lassosum2.sscore"))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_tune <- read.delim("All_Tune.txt")
file.remove("All_Tune.txt")
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[,-c(1:27)]

#Validation
prs_validation_CT <- read.csv(paste0(trait,"_prs_all_validation.txt"), sep="")
prs_validation_LDPred2 <- read.delim(paste0(trait,"_prs_validation_ldpred2.sscore"))
prs_validation_LASSOSum2 <- read.delim(paste0(trait,"_prs_validation_lassosum2.sscore"))

file.remove(paste0(trait,"_prs_all_validation.txt"))
file.remove(paste0(trait,"_prs_validation_ldpred2.sscore"))
file.remove(paste0(trait,"_prs_validation_lassosum2.sscore"))

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],prs_validation_LDPred2[,-c(1,2,3,4)],prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_validation <- read.delim("All_Validation.txt")
file.remove("All_Validation.txt")
pheno_validation <- left_join(pheno_validation,prs_validation_all,by = "IID")

prs_validation_all <- pheno_validation[,-c(1:27)]

## Null Models

model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_train)
y_train <- model.null$residual

model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
y_tune <- model.null$residual

model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
y_validation <- model.null$residual

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

pheno_tune$y_tune <- NA
pheno_tune$y_tune[!is.na(pheno_tune[,trait])] <- y_tune

pheno_validation$y_validation <- NA
pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- y_validation

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

Results <- Ensemble_Function(x = prs_tune_all,y = pheno_tune[,"y_tune"],family = "continuous")
Results$Coefficients[is.na(Results$Coefficients)] <- 0
write.csv(Results$Coefficients,file = paste0(trait,"_Final_Coefficients.csv"))

PRS_Train <- as.matrix(pheno_train[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Tune <- as.matrix(pheno_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Validation <- as.matrix(pheno_validation[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)

prs_best_train <- data.frame(IID = pheno_train$IID,prs = PRS_Train)
prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = PRS_Tune)
prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = PRS_Validation)

write.table(prs_best_train,file=paste0(trait,"_Best_Train_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0(trait,"_Best_Tune_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0(trait,"_Best_Validation_All.txt"),sep = "\t",row.names = FALSE)

load("all_phenotypes.RData")
file.remove("all_phenotypes.RData")


pheno_validation_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),"y_validation")],prs_best_validation)
pheno_validation_adjusted <- pheno_validation_raw


tmp <- data.frame(y = pheno_validation_adjusted[,"prs"],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
R <- mod$residuals
tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
y_hat <- predict(mod,tmp)
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

pheno_validation_raw_EUR[,"prs"] <- scale(pheno_validation_raw_EUR[,"prs"])
pheno_validation_raw_SAS[,"prs"] <- scale(pheno_validation_raw_SAS[,"prs"])
pheno_validation_raw_AMR[,"prs"] <- scale(pheno_validation_raw_AMR[,"prs"])
pheno_validation_raw_AFR[,"prs"] <- scale(pheno_validation_raw_AFR[,"prs"])
pheno_validation_raw_EAS[,"prs"] <- scale(pheno_validation_raw_EAS[,"prs"])

Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(as.formula(paste0("y_validation~","prs")),data = boot_data))[2]
  return(c(result))
}

R2_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- summary(lm(as.formula(paste0("y_validation~","prs")),data = boot_data))$r.squared
  return(c(result))
}

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_EUR))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 10000)
beta_raw_EUR_boot <- boot_beta$t
beta_se_validation_raw_EUR <- sd(boot_beta$t)

R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 10000)
R2_raw_EUR_boot <- boot_R2$t
R2_se_validation_raw_EUR <- sd(boot_R2$t)

beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_SAS))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 10000)
beta_raw_SAS_boot <- boot_beta$t
beta_se_validation_raw_SAS <- sd(boot_beta$t)

R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 10000)
R2_raw_SAS_boot <- boot_R2$t
R2_se_validation_raw_SAS <- sd(boot_R2$t)

beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_AMR))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 10000)
beta_raw_AMR_boot <- boot_beta$t
beta_se_validation_raw_AMR <- sd(boot_beta$t)

R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 10000)
R2_raw_AMR_boot <- boot_R2$t
R2_se_validation_raw_AMR <- sd(boot_R2$t)

beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_AFR))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 10000)
beta_raw_AFR_boot <- boot_beta$t
beta_se_validation_raw_AFR <- sd(boot_beta$t)

R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 10000)
R2_raw_AFR_boot <- boot_R2$t
R2_se_validation_raw_AFR <- sd(boot_R2$t)

beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_EAS))[2]
boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 10000)
beta_raw_EAS_boot <- boot_beta$t
beta_se_validation_raw_EAS <- sd(boot_beta$t)

R2_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_raw_EAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EAS, statistic = R2_Boot, R = 10000)
R2_raw_EAS_boot <- boot_R2$t
R2_se_validation_raw_EAS <- sd(boot_R2$t)

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_EUR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 10000)
beta_adjusted_EUR_boot <- boot_beta$t
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)

R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 10000)
R2_adjusted_EUR_boot <- boot_R2$t
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)

beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_SAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 10000)
beta_adjusted_SAS_boot <- boot_beta$t
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)

R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 10000)
R2_adjusted_SAS_boot <- boot_R2$t
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)

beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_AMR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 10000)
beta_adjusted_AMR_boot <- boot_beta$t
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)

R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 10000)
R2_adjusted_AMR_boot <- boot_R2$t
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)

beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_AFR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 10000)
beta_adjusted_AFR_boot <- boot_beta$t
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)

R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 10000)
R2_adjusted_AFR_boot <- boot_R2$t
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)

beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_EAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 10000)
beta_adjusted_EAS_boot <- boot_beta$t
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)

R2_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~","prs")),data = pheno_validation_adjusted_EAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EAS, statistic = R2_Boot, R = 10000)
R2_adjusted_EAS_boot <- boot_R2$t
R2_se_validation_adjusted_EAS <- sd(boot_R2$t)

CV_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                         beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                         R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                         R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                         beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                         beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                         R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                         R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))

CV_Boot_Results <- data.frame(trait = trait,beta_raw_EUR_boot,R2_raw_EUR_boot,beta_raw_SAS_boot,R2_raw_SAS_boot,
                              beta_raw_AMR_boot,R2_raw_AMR_boot,beta_raw_AFR_boot,R2_raw_AFR_boot,
                              beta_raw_EAS_boot,R2_raw_EAS_boot,beta_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                              beta_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                              beta_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_adjusted_EAS_boot,R2_adjusted_EAS_boot)

write.csv(CV_Results,file = paste0(trait,"Best_Betas.csv"),row.names = FALSE) 
write.csv(CV_Boot_Results,file = paste0(trait,"_Bootstraps_Results.csv"),row.names = FALSE)  