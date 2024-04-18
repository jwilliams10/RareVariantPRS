rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(stringr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

#Train
prs_train_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_train",i,".txt"), sep="")
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_train",i,".txt")))
prs_train_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/prs_train",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/prs_train",i,".sscore")))
prs_train_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/prs_train",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/prs_train",i,".sscore")))

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],prs_train_LDPred2[,-c(1,2,3,4)],prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")
## Merge covariates and y for tuning with the prs_mat
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[,-c(1:3)]

#Tune
prs_tune_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_tune",i,".txt"), sep="")
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_tune",i,".txt")))
prs_tune_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/prs_tune",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/prs_tune",i,".sscore")))
prs_tune_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/prs_tune",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/prs_tune",i,".sscore")))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[,-c(1:3)]

#Validation
prs_validation_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_validation",i,".txt"), sep="")
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_validation",i,".txt")))
prs_validation_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/prs_validation",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/prs_validation",i,".sscore")))
prs_validation_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/prs_validation",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/prs_validation",i,".sscore")))

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],prs_validation_LDPred2[,-c(1,2,3,4)],prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")
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

## SL

SL.library <- c(
  "SL.glmnet",
  "SL.ridge",
  "SL.glm",
  "SL.mean"
)
sl <- SuperLearner(Y = y_tune, X = prs_tune_all, family = gaussian(),
                   # For a real analysis we would use V = 10.
                   # V = 3,
                   SL.library = SL.library,cvControl = list(V = 20))
cvsl <- CV.SuperLearner(Y = y_tune, X = prs_tune_all, family = gaussian(),
                        # For a real analysis we would use V = 10.
                        # V = 3,
                        SL.library = SL.library,cvControl = list(V = 20))

best_algorithm <- summary(cvsl)$Table$Algorithm[which.min(summary(cvsl)$Table$Ave)]

a <- predict(sl, prs_validation_all, onlySL = FALSE)

prs_best_validation_sl <- a$pred
prs_best_validation_glmnet <- a$library.predict[,1]
prs_best_validation_ridge <- a$library.predict[,2]
prs_best_validation_glm <- a$library.predict[,3]

if(best_algorithm == "SL.glmnet_All"){
  #final
  prs_best_validation <- prs_best_validation_glmnet
}else if(best_algorithm == "SL.ridge_All"){
  #final
  prs_best_validation <- prs_best_validation_ridge
}else if(best_algorithm == "SL.glm_All"){
  #final
  prs_best_validation <- prs_best_validation_glm
}else{
  #final
  prs_best_validation <- prs_best_validation_sl
}

prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = prs_best_validation)


a <- predict(sl, prs_train_all, onlySL = FALSE)

prs_best_train_sl <- a$pred
prs_best_train_glmnet <- a$library.predict[,1]
prs_best_train_ridge <- a$library.predict[,2]
prs_best_train_glm <- a$library.predict[,3]

if(best_algorithm == "SL.glmnet_All"){
  #final
  prs_best_train <- prs_best_train_glmnet
}else if(best_algorithm == "SL.ridge_All"){
  #final
  prs_best_train <- prs_best_train_ridge
}else if(best_algorithm == "SL.glm_All"){
  #final
  prs_best_train <- prs_best_train_glm
}else{
  #final
  prs_best_train <- prs_best_train_sl
}
prs_best_train <- data.frame(IID = pheno_train$IID,prs = prs_best_train)


a <- predict(sl, prs_tune_all, onlySL = FALSE)

prs_best_tune_sl <- a$pred
prs_best_tune_glmnet <- a$library.predict[,1]
prs_best_tune_ridge <- a$library.predict[,2]
prs_best_tune_glm <- a$library.predict[,3]

if(best_algorithm == "SL.glmnet_All"){
  #final
  prs_best_tune <- prs_best_tune_glmnet
}else if(best_algorithm == "SL.ridge_All"){
  #final
  prs_best_tune <- prs_best_tune_ridge
}else if(best_algorithm == "SL.glm_All"){
  #final
  prs_best_tune <- prs_best_tune_glm
}else{
  #final
  prs_best_tune <- prs_best_tune_sl
}
prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = prs_best_tune)

pheno_train <- left_join(pheno_train,prs_best_train)
pheno_validation <- left_join(pheno_validation,prs_best_validation)
pheno_tune <- left_join(pheno_tune,prs_best_tune)

prs_columns <- c(which(str_detect(colnames(pheno_tune),"CT_")),which(str_detect(colnames(pheno_tune),"LDPred2_")),which(str_detect(colnames(pheno_tune),"LASSOSum2_")),which(str_detect(colnames(pheno_tune),"prs")))

r2_tune <- vector()
for(j in 1:length(prs_columns)){
  r2_tune[j] <- summary(lm(as.formula(paste0("y_tune ~",colnames(pheno_tune)[prs_columns[j]])),data = pheno_tune))$r.squared
}



prs_best_train <- data.frame(IID = pheno_train$IID,prs = pheno_train[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = pheno_tune[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = pheno_validation[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

write.table(prs_best_train,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Train_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")


pheno_validation_raw <- pheno_validation

pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_raw_NonEUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_raw_UNK <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_MIX <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
pheno_validation_raw_NonEUR$y_validation <- scale(pheno_validation_raw_NonEUR$y_validation)
pheno_validation_raw_UNK$y_validation <- scale(pheno_validation_raw_UNK$y_validation)
pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
pheno_validation_raw_MIX$y_validation <- scale(pheno_validation_raw_MIX$y_validation)
pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)
pheno_validation_raw_EAS$y_validation <- scale(pheno_validation_raw_EAS$y_validation)

pheno_validation_raw_EUR[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_EUR[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])
pheno_validation_raw_NonEUR[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_NonEUR[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])
pheno_validation_raw_UNK[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_UNK[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])
pheno_validation_raw_SAS[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_SAS[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])
pheno_validation_raw_MIX[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_MIX[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])
pheno_validation_raw_AFR[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_AFR[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])
pheno_validation_raw_EAS[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- scale(pheno_validation_raw_EAS[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_EUR))[2]
se_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_EUR))$coefficients[2,2]
beta_validation_raw_NonEUR <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_NonEUR))[2]
se_validation_raw_NonEUR <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_NonEUR))$coefficients[2,2]
beta_validation_raw_UNK <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_UNK))[2]
se_validation_raw_UNK <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_UNK))$coefficients[2,2]
beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_SAS))[2]
se_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_SAS))$coefficients[2,2]
beta_validation_raw_MIX <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_MIX))[2]
se_validation_raw_MIX <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_MIX))$coefficients[2,2]
beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_AFR))[2]
se_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_AFR))$coefficients[2,2]
beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_EAS))[2]
se_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_raw_EAS))$coefficients[2,2]

CV_Results <- data.frame(i = i,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_NonEUR,beta_validation_raw_UNK,beta_validation_raw_SAS,beta_validation_raw_MIX,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                         se_raw = c(se_validation_raw_EUR,se_validation_raw_NonEUR,se_validation_raw_UNK,se_validation_raw_SAS,se_validation_raw_MIX,se_validation_raw_AFR,se_validation_raw_EAS))

write.csv(CV_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Betas",i,".csv"),row.names = FALSE)