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

prs_train_all <- pheno_train[!is.na(pheno_train[,trait]),-c(1:27)]

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

prs_tune_all <- pheno_tune[!is.na(pheno_tune[,trait]),-c(1:27)]

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

prs_validation_all <- pheno_validation[!is.na(pheno_validation[,trait]),-c(1:27)]

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

prs_best_validation <- data.frame(IID = pheno_validation$IID[!is.na(pheno_validation[,trait])],prs = prs_best_validation)


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
prs_best_train <- data.frame(IID = pheno_train$IID[!is.na(pheno_train[,trait])],prs = prs_best_train)


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
prs_best_tune <- data.frame(IID = pheno_tune$IID[!is.na(pheno_tune[,trait])],prs = prs_best_tune)

pheno_train <- left_join(pheno_train,prs_best_train)
pheno_validation <- left_join(pheno_validation,prs_best_validation)
pheno_tune <- left_join(pheno_tune,prs_best_tune)

prs_columns <- c(which(str_detect(colnames(pheno_tune),"CT_")),which(str_detect(colnames(pheno_tune),"LDPred2_")),which(str_detect(colnames(pheno_tune),"LASSOSum2_")),which(str_detect(colnames(pheno_tune),"prs")))

r2_tune <- vector()
for(i in 1:length(prs_columns)){
  r2_tune[i] <- summary(lm(as.formula(paste0("y_tune ~",colnames(pheno_tune)[prs_columns[i]])),data = pheno_tune))$r.squared
}

if(colnames(pheno_tune)[prs_columns[which.max(r2_tune)]] == "prs"){
  if(best_algorithm == "SL.glmnet_All"){
    beta_sl <- data.frame(Coef = row.names(coef(sl$fitLibrary$SL.glmnet_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glmnet_All$object)))
  }else if(best_algorithm == "SL.ridge_All"){
    beta_sl <- data.frame(Coef = row.names(sl$fitLibrary$SL.ridge_All$bestCoef),Beta = as.numeric(sl$fitLibrary$SL.ridge_All$bestCoef[,1]))
  }else if(best_algorithm == "SL.glm_All"){
    beta_sl <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glm_All$object)))
  }else{
    beta_lasso <- data.frame(Coef = row.names(coef(sl$fitLibrary$SL.glmnet_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glmnet_All$object)))
    beta_ridge <- data.frame(Coef = row.names(sl$fitLibrary$SL.ridge_All$bestCoef),Beta = as.numeric(sl$fitLibrary$SL.ridge_All$bestCoef[,1]))
    beta_glm <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)),Beta = as.numeric(coef(sl$fitLibrary$SL.glm_All$object)))
    beta_mean <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)), Beta = 0)
    beta_mean$Beta[1] <- sl$fitLibrary$SL.mean_All$object
    
    beta_sl <- data.frame(Coef = names(coef(sl$fitLibrary$SL.glm_All$object)), Beta = 0)
    beta_sl$Beta <- as.numeric(sl$coef[names(sl$coef) == "SL.glmnet_All"]) * beta_lasso$Beta + 
      as.numeric(sl$coef[names(sl$coef) == "SL.ridge_All"]) * beta_ridge$Beta + 
      as.numeric(sl$coef[names(sl$coef) == "SL.glm_All"]) * beta_glm$Beta + 
      as.numeric(sl$coef[names(sl$coef) == "SL.mean_All"]) * beta_mean$Beta 
  }
  write.csv(beta_sl,file = paste0(trait,"_Final_Coefficients.csv"),row.names = FALSE)
  
}else{
  beta_coefs <- data.frame(Coef = colnames(pheno_tune)[prs_columns[which.max(r2_tune)]],Beta = 1)
  write.csv(beta_coefs,file = paste0(trait,"_Final_Coefficients.csv"),row.names = FALSE)
}






prs_best_train <- data.frame(IID = pheno_train$IID,prs = pheno_train[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = pheno_tune[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = pheno_validation[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

write.table(prs_best_train,file=paste0(trait,"_Best_Train_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0(trait,"_Best_Tune_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0(trait,"_Best_Validation_All.txt"),sep = "\t",row.names = FALSE)

load("all_phenotypes.RData")
file.remove("all_phenotypes.RData")


pheno_validation_raw <- pheno_validation
pheno_validation_adjusted <- pheno_validation


tmp <- data.frame(y = pheno_validation_adjusted[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]],pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
R <- mod$residuals
tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
y_hat <- predict(mod,tmp)
if(sum(sqrt(y_hat)) == 0){
  pheno_validation_adjusted[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- 0
}else{
  pheno_validation_adjusted[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]] <- R/sqrt(y_hat)
}


pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_raw_NonEUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_raw_UNK <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_MIX <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_NonEUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_adjusted_UNK <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_MIX <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
pheno_validation_raw_NonEUR$y_validation <- scale(pheno_validation_raw_NonEUR$y_validation)
pheno_validation_raw_UNK$y_validation <- scale(pheno_validation_raw_UNK$y_validation)
pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
pheno_validation_raw_MIX$y_validation <- scale(pheno_validation_raw_MIX$y_validation)
pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)
pheno_validation_raw_EAS$y_validation <- scale(pheno_validation_raw_EAS$y_validation)

pheno_validation_adjusted_EUR$y_validation <- scale(pheno_validation_adjusted_EUR$y_validation)
pheno_validation_adjusted_NonEUR$y_validation <- scale(pheno_validation_adjusted_NonEUR$y_validation)
pheno_validation_adjusted_UNK$y_validation <- scale(pheno_validation_adjusted_UNK$y_validation)
pheno_validation_adjusted_SAS$y_validation <- scale(pheno_validation_adjusted_SAS$y_validation)
pheno_validation_adjusted_MIX$y_validation <- scale(pheno_validation_adjusted_MIX$y_validation)
pheno_validation_adjusted_AFR$y_validation <- scale(pheno_validation_adjusted_AFR$y_validation)
pheno_validation_adjusted_EAS$y_validation <- scale(pheno_validation_adjusted_EAS$y_validation)


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

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_EUR))[2]
se_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_EUR))$coefficients[2,2]
beta_validation_adjusted_NonEUR <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_NonEUR))[2]
se_validation_adjusted_NonEUR <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_NonEUR))$coefficients[2,2]
beta_validation_adjusted_UNK <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_UNK))[2]
se_validation_adjusted_UNK <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_UNK))$coefficients[2,2]
beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_SAS))[2]
se_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_SAS))$coefficients[2,2]
beta_validation_adjusted_MIX <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_MIX))[2]
se_validation_adjusted_MIX <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_MIX))$coefficients[2,2]
beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_AFR))[2]
se_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_AFR))$coefficients[2,2]
beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_EAS))[2]
se_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~",colnames(pheno_tune)[prs_columns[which.max(r2_tune)]])),data = pheno_validation_adjusted_EAS))$coefficients[2,2]

CV_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_NonEUR,beta_validation_raw_UNK,beta_validation_raw_SAS,beta_validation_raw_MIX,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                         se_raw = c(se_validation_raw_EUR,se_validation_raw_NonEUR,se_validation_raw_UNK,se_validation_raw_SAS,se_validation_raw_MIX,se_validation_raw_AFR,se_validation_raw_EAS), 
                         beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_NonEUR,beta_validation_adjusted_UNK,beta_validation_adjusted_SAS,beta_validation_adjusted_MIX,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                         se_adjusted = c(se_validation_adjusted_EUR,se_validation_adjusted_NonEUR,se_validation_adjusted_UNK,se_validation_adjusted_SAS,se_validation_adjusted_MIX,se_validation_adjusted_AFR,se_validation_adjusted_EAS))

write.csv(CV_Results,file = paste0(trait,"Best_Betas.csv"),row.names = FALSE)