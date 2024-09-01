rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)
library(stringr)

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
prs_train_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_train.txt"), sep="")
prs_train_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_train.sscore"))
prs_train_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_train.sscore"))

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],(-1)*prs_train_LDPred2[,-c(1,2,3,4)],(-1)*prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[!is.na(pheno_train[,trait]),-c(1:27)]
pheno_train <- pheno_train[!is.na(pheno_train[,trait]),]

#Tune
prs_tune_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
prs_tune_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_tune.sscore"))
prs_tune_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_tune.sscore"))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],(-1)*prs_tune_LDPred2[,-c(1,2,3,4)],(-1)*prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[!is.na(pheno_tune[,trait]),-c(1:27)]
pheno_tune <- pheno_tune[!is.na(pheno_tune[,trait]),]

#Validation
prs_validation_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
prs_validation_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_validation.sscore"))
prs_validation_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_validation.sscore"))

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],(-1)*prs_validation_LDPred2[,-c(1,2,3,4)],(-1)*prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
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


## SL

SL.library <- c(
  "SL.glmnet",
  "SL.glm",
  "SL.mean"
)
sl <- SuperLearner(Y = pheno_tune[,trait], X = prs_tune_all, family = binomial(), method = "method.AUC",
                   # For a real analysis we would use V = 10.
                   # V = 3,
                   SL.library = SL.library,cvControl = list(V = 20,stratifyCV = TRUE))
cvsl <- CV.SuperLearner(Y = pheno_tune[,trait], X = prs_tune_all, family = binomial(), method = "method.AUC",
                        # For a real analysis we would use V = 10.
                        # V = 3,
                        SL.library = SL.library,cvControl = list(V = 20,stratifyCV = TRUE))

best_algorithm <- summary(cvsl)$Table$Algorithm[which.max(summary(cvsl)$Table$Ave)]

a <- predict(sl, prs_validation_all, onlySL = FALSE)

prs_best_validation_sl <- log(a$pred)/(1 - log(a$pred))
prs_best_validation_glmnet <- log(a$library.predict[,1])/(1 - log(a$library.predict[,1]))
prs_best_validation_ridge <- log(a$library.predict[,2])/(1 - log(a$library.predict[,2]))
prs_best_validation_glm <- log(a$library.predict[,3])/(1 - log(a$library.predict[,3]))

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

prs_best_train_sl <- log(a$pred)/(1 - log(a$pred))
prs_best_train_glmnet <- log(a$library.predict[,1])/(1 - log(a$library.predict[,1]))
prs_best_train_ridge <- log(a$library.predict[,2])/(1 - log(a$library.predict[,2]))
prs_best_train_glm <- log(a$library.predict[,3])/(1 - log(a$library.predict[,3]))

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
prs_best_train <- data.frame(IID = pheno_train$IID[!is.na(pheno_train[,trait])],prs = prs_best_train)


a <- predict(sl, prs_tune_all, onlySL = FALSE)

prs_best_tune_sl <- log(a$pred)/(1 - log(a$pred))
prs_best_tune_glmnet <- log(a$library.predict[,1])/(1 - log(a$library.predict[,1]))
prs_best_tune_ridge <- log(a$library.predict[,2])/(1 - log(a$library.predict[,2]))
prs_best_tune_glm <- log(a$library.predict[,3])/(1 - log(a$library.predict[,3]))

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
prs_best_tune <- data.frame(IID = pheno_tune$IID[!is.na(pheno_tune[,trait])],prs = prs_best_tune)

pheno_train <- left_join(pheno_train,prs_best_train)
pheno_validation <- left_join(pheno_validation,prs_best_validation)
pheno_tune <- left_join(pheno_tune,prs_best_tune)

if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <- paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}

prs_columns <- c(which(str_detect(colnames(pheno_tune),"CT_")),which(str_detect(colnames(pheno_tune),"LDPred2_")),which(str_detect(colnames(pheno_tune),"LASSOSum2_")),which(str_detect(colnames(pheno_tune),"prs")))

auc_tune <- vector()
for(i in 1:length(prs_columns)){
  
  auc_tune[i] <- roc.binary(status = trait,
                            variable = colnames(pheno_tune)[prs_columns[i]],
                            confounders = as.formula(confounders),
                            data = pheno_tune[!is.na(pheno_tune[,trait]),],
                            precision=seq(0.05,0.95, by=0.05))$auc
}


if(colnames(pheno_tune)[prs_columns[which.max(auc_tune)]] == "prs"){
  if(best_algorithm == "SL.glmnet_All"){
    beta_sl <- data.frame(Coef = row.names(coef(sl$fitLibrary$SL.glmnet_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glmnet_All$object)))
  }else if(best_algorithm == "SL.glm_All"){
    beta_sl <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glm_All$object)))
  }else{
    beta_lasso <- data.frame(Coef = row.names(coef(sl$fitLibrary$SL.glmnet_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glmnet_All$object)))
    beta_glm <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glm_All$object)))
    beta_mean <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)), Beta = 0)
    beta_mean$Beta[1] <- sl$fitLibrary$SL.mean_All$object
    
    beta_sl <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)), Beta = 0)
    beta_sl$Beta <- as.numeric(sl$coef[names(sl$coef) == "SL.glmnet_All"]) * beta_lasso$Beta + 
      as.numeric(sl$coef[names(sl$coef) == "SL.glm_All"]) * beta_glm$Beta + 
      as.numeric(sl$coef[names(sl$coef) == "SL.mean_All"]) * beta_mean$Beta 
  }
  write.csv(beta_sl,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_final_coef.csv"),row.names = FALSE)
  
}else{
  beta_coefs <- data.frame(Coef = colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],Beta = 1)
  write.csv(beta_coefs,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_final_coef.csv"),row.names = FALSE)
}






prs_best_train <- data.frame(IID = pheno_train$IID,prs = pheno_train[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_train <- pheno_train[,colnames(pheno_train) != "prs"]
pheno_train <- left_join(pheno_train,prs_best_train)

prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = pheno_tune[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_tune <- pheno_tune[,colnames(pheno_tune) != "prs"]
pheno_tune <- left_join(pheno_tune,prs_best_tune)

prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = pheno_validation[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_validation <- pheno_validation[,colnames(pheno_validation) != "prs"]
pheno_validation <- left_join(pheno_validation,prs_best_validation)

write.table(prs_best_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Train_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"),sep = "\t",row.names = FALSE)




load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")


pheno_validation_raw <- pheno_validation
pheno_validation_adjusted <- pheno_validation


tmp <- data.frame(y = pheno_validation_adjusted[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
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
  pheno_validation_adjusted[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- 0
}else{
  pheno_validation_adjusted[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- R/sqrt(y_hat)
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

pheno_validation_raw_EUR[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- scale(pheno_validation_raw_EUR[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_validation_raw_SAS[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- scale(pheno_validation_raw_SAS[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_validation_raw_AMR[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- scale(pheno_validation_raw_AMR[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_validation_raw_AFR[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- scale(pheno_validation_raw_AFR[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])
pheno_validation_raw_EAS[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]] <- scale(pheno_validation_raw_EAS[,colnames(pheno_tune)[prs_columns[which.max(auc_tune)]]])


beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
se_validation_raw_EUR <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))$coefficients[2,2]
beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
se_validation_raw_SAS <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))$coefficients[2,2]
beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
se_validation_raw_AMR <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))$coefficients[2,2]
beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
se_validation_raw_AFR <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))$coefficients[2,2]
beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
se_validation_raw_EAS <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))$coefficients[2,2]

beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
se_validation_adjusted_EUR <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
se_validation_adjusted_SAS <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
se_validation_adjusted_AMR <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
se_validation_adjusted_AFR <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
se_validation_adjusted_EAS <- summary(glm(as.formula(paste0(trait,"~",colnames(pheno_tune)[prs_columns[which.max(auc_tune)]],"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))$coefficients[2,2]

SL_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                         se_raw = c(se_validation_raw_EUR,se_validation_raw_SAS,se_validation_raw_AMR,se_validation_raw_AFR,se_validation_raw_EAS), 
                         beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                         se_adjusted = c(se_validation_adjusted_EUR,se_validation_adjusted_SAS,se_validation_adjusted_AMR,se_validation_adjusted_AFR,se_validation_adjusted_EAS))

write.csv(SL_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"Best_Betas.csv"),row.names = FALSE)