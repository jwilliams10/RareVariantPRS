rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(stringr)
library(glmnet)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

#Train
prs_train_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/prs_all_train",i,".txt"), sep="")
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/prs_all_train",i,".txt")))
prs_train_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/prs_train",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/prs_train",i,".sscore")))
prs_train_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/prs_train",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/prs_train",i,".sscore")))

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],prs_train_LDPred2[,-c(1,2,3,4)],prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")
## Merge covariates and y for tuning with the prs_mat
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[,-c(1:3)]

#Tune
prs_tune_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/prs_all_tune",i,".txt"), sep="")
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/prs_all_tune",i,".txt")))
prs_tune_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/prs_tune",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/prs_tune",i,".sscore")))
prs_tune_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/prs_tune",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/prs_tune",i,".sscore")))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[,-c(1:3)]

#Validation
prs_validation_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/prs_all_validation",i,".txt"), sep="")
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/prs_all_validation",i,".txt")))
prs_validation_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/prs_validation",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/prs_validation",i,".sscore")))
prs_validation_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/prs_validation",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/prs_validation",i,".sscore")))

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],prs_validation_LDPred2[,-c(1,2,3,4)],prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/phenotypes/Y_Validation.RData")
pheno_validation <- Y_validation[[i]]
colnames(pheno_validation) <- c("IID","Y")
pheno_validation <- left_join(pheno_validation,prs_validation_all,by = "IID")

prs_validation_all <- pheno_validation[,-c(1:3)]

## Null Models

model.null <- lm(Y~1,data=pheno_train)
y_train <- model.null$residual

model.null <- lm(Y~1,data=pheno_tune)
y_tune <- model.null$residual

model.null <- lm(Y~1,data=pheno_validation)
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
pheno_tune$y_tune <- y_tune

pheno_validation$y_validation <- NA
pheno_validation$y_validation <- y_validation


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

## SL

Results <- Ensemble_Function(x = prs_tune_all,y = pheno_tune[,"y_tune"],family = "continuous")
Results$Coefficients[is.na(Results$Coefficients)] <- 0

PRS_Train <- as.matrix(pheno_train[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Tune <- as.matrix(pheno_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Validation <- as.matrix(pheno_validation[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)

PRS_Train <- data.frame(IID = pheno_train$IID,PRS = PRS_Train)
PRS_Tune <- data.frame(IID = pheno_tune$IID,PRS = PRS_Tune)
PRS_Validation <- data.frame(IID = pheno_validation$IID,PRS = PRS_Validation)

write.table(PRS_Train,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Combined_Common_PRS/Best_Train_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(PRS_Tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(PRS_Validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")


pheno_validation_raw <- inner_join(pheno_validation[,c("IID","y_validation")],PRS_Validation)
pheno_validation_adjusted <- pheno_validation_raw
pheno_validation_adjusted <- inner_join(pheno_validation_adjusted,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])

tmp <- data.frame(y = pheno_validation_adjusted[,"PRS"],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
R <- mod$residuals
tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
y_hat <- predict(mod,tmp)
if(sum(sqrt(y_hat)) == 0){
  pheno_validation_adjusted[,"PRS"] <- 0
}else{
  pheno_validation_adjusted[,"PRS"] <- R/sqrt(y_hat)
}


pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]

pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
pheno_validation_raw_AMR$y_validation <- scale(pheno_validation_raw_AMR$y_validation)
pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)

pheno_validation_adjusted_EUR$y_validation <- scale(pheno_validation_adjusted_EUR$y_validation)
pheno_validation_adjusted_SAS$y_validation <- scale(pheno_validation_adjusted_SAS$y_validation)
pheno_validation_adjusted_AMR$y_validation <- scale(pheno_validation_adjusted_AMR$y_validation)
pheno_validation_adjusted_AFR$y_validation <- scale(pheno_validation_adjusted_AFR$y_validation)


pheno_validation_raw_EUR[,"PRS"] <- scale(pheno_validation_raw_EUR[,"PRS"])
pheno_validation_raw_SAS[,"PRS"] <- scale(pheno_validation_raw_SAS[,"PRS"])
pheno_validation_raw_AMR[,"PRS"] <- scale(pheno_validation_raw_AMR[,"PRS"])
pheno_validation_raw_AFR[,"PRS"] <- scale(pheno_validation_raw_AFR[,"PRS"])

Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = boot_data))[2]
  return(c(result))
}

R2_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = boot_data))$r.squared
  return(c(result))
}

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_EUR))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EUR <- sd(boot_R2$t)
R2_lower_validation_raw_EUR <- R2_ci$basic[4]
R2_upper_validation_raw_EUR <- R2_ci$basic[5]

beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_SAS))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_SAS <- sd(boot_R2$t)
R2_lower_validation_raw_SAS <- R2_ci$basic[4]
R2_upper_validation_raw_SAS <- R2_ci$basic[5]

beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AMR))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AMR <- sd(boot_R2$t)
R2_lower_validation_raw_AMR <- R2_ci$basic[4]
R2_upper_validation_raw_AMR <- R2_ci$basic[5]

beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AFR))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_raw_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AFR <- sd(boot_R2$t)
R2_lower_validation_raw_AFR <- R2_ci$basic[4]
R2_upper_validation_raw_AFR <- R2_ci$basic[5]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_EUR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
R2_lower_validation_adjusted_EUR <- R2_ci$basic[4]
R2_upper_validation_adjusted_EUR <- R2_ci$basic[5]

beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_SAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_SAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_SAS <- R2_ci$basic[5]

beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AMR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AMR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AMR <- R2_ci$basic[5]

beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AFR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = pheno_validation_adjusted_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AFR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AFR <- R2_ci$basic[5]

RICE_CV <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR"), 
                      beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR), 
                      beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR), 
                      beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR), 
                      beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR), 
                      R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR),
                      R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR),
                      R2_lower_raw = c(R2_lower_validation_raw_EUR,R2_lower_validation_raw_SAS,R2_lower_validation_raw_AMR,R2_lower_validation_raw_AFR),
                      R2_upper_raw = c(R2_upper_validation_raw_EUR,R2_upper_validation_raw_SAS,R2_upper_validation_raw_AMR,R2_upper_validation_raw_AFR),
                      beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR), 
                      beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR), 
                      beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR), 
                      beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR), 
                      R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR),
                      R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR),
                      R2_lower_adjusted = c(R2_lower_validation_adjusted_EUR,R2_lower_validation_adjusted_SAS,R2_lower_validation_adjusted_AMR,R2_lower_validation_adjusted_AFR),
                      R2_upper_adjusted = c(R2_upper_validation_adjusted_EUR,R2_upper_validation_adjusted_SAS,R2_upper_validation_adjusted_AMR,R2_upper_validation_adjusted_AFR))


write.csv(RICE_CV,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Combined_Common_PRS/Best_Betas",i,".csv"),row.names = FALSE)