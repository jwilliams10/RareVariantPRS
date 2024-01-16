rm(list = ls())
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
prs_train_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_train.txt"), sep="")
prs_train_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_prs_train.sscore"))
prs_train_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_prs_train.sscore"))

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],prs_train_LDPred2[,-c(1,2,3,4)],prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[!is.na(pheno_train[,trait]),-c(1:27)]

#Tune
prs_tune_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
prs_tune_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_prs_tune.sscore"))
prs_tune_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_prs_tune.sscore"))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[!is.na(pheno_tune[,trait]),-c(1:27)]

#Validation
prs_validation_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
prs_validation_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_prs_validation.sscore"))
prs_validation_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_prs_validation.sscore"))

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],prs_validation_LDPred2[,-c(1,2,3,4)],prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
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
prs_best_train <- data.frame(IID = pheno_train$IID,prs = pheno_train[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = pheno_tune[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = pheno_validation[,colnames(pheno_tune)[prs_columns[which.max(r2_tune)]]])

write.table(prs_best_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Train_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

y_validation_EUR <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]]
y_validation_NonEUR <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"]]
y_validation_UNK <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"]]
y_validation_SAS <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]]
y_validation_MIX <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"]]
y_validation_AFR <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]]
y_validation_EAS <- y_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]]

prs_best_validation_EUR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
prs_best_validation_NonEur <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
prs_best_validation_UNK <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
prs_best_validation_SAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
prs_best_validation_MIX <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
prs_best_validation_AFR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
prs_best_validation_EAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_EUR~prs_best_validation_EUR$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_EUR, x = prs_best_validation_EUR$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_Eur",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur.RData"))

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_NonEUR~prs_best_validation_NonEur$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_NonEUR, x = prs_best_validation_NonEur$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_NonEur",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur.RData"))

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_UNK~prs_best_validation_UNK$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_UNK, x = prs_best_validation_UNK$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK.RData"))

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_SAS~prs_best_validation_SAS$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_SAS, x = prs_best_validation_SAS$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_SAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS.RData"))

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_MIX~prs_best_validation_MIX$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_MIX, x = prs_best_validation_MIX$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_MIX",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX.RData"))

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_AFR~prs_best_validation_AFR$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_AFR, x = prs_best_validation_AFR$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_AFR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR.RData"))

## bootstrap the R2 to provide an approximate distribution 
model <- lm(y_validation_EAS~prs_best_validation_EAS$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_validation_EAS, x = prs_best_validation_EAS$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined_EAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS.RData"))

