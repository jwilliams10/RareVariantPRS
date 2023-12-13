rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

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

prs_validation_all <- pheno_validation[!is.na(pheno_validation[,trait]),-c(1:27)]
pheno_validation <- pheno_validation[!is.na(pheno_validation[,trait]),]

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
#  "SL.ridge",
  "SL.glm",
  "SL.mean"
)
sl <- SuperLearner(Y = pheno_tune[,trait], X = prs_tune_all, family = binomial(), method = "method.AUC",
                   # For a real analysis we would use V = 10.
                   # V = 3,
                   SL.library = SL.library)
cvsl <- CV.SuperLearner(Y = pheno_tune[,trait], X = prs_tune_all, family = binomial(), method = "method.AUC",
                        # For a real analysis we would use V = 10.
                        # V = 3,
                        SL.library = SL.library)

best_algorithm <- summary(cvsl)$Table$Algorithm[which.max(summary(cvsl)$Table$Ave)]

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

write.table(prs_best_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Train_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad_EUR <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEUR <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

prs_best_validation_EUR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
prs_best_validation_NonEur <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
prs_best_validation_UNK <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
prs_best_validation_SAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
prs_best_validation_MIX <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
prs_best_validation_AFR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
prs_best_validation_EAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_vad_EUR <- inner_join(pheno_vad_EUR,prs_best_validation_EUR)
pheno_vad_NonEUR <- inner_join(pheno_vad_NonEUR,prs_best_validation_NonEur)
pheno_vad_UNK <- inner_join(pheno_vad_UNK,prs_best_validation_UNK)
pheno_vad_SAS <- inner_join(pheno_vad_SAS,prs_best_validation_SAS)
pheno_vad_MIX <- inner_join(pheno_vad_MIX,prs_best_validation_MIX)
pheno_vad_AFR <- inner_join(pheno_vad_AFR,prs_best_validation_AFR)
pheno_vad_EAS <- inner_join(pheno_vad_EAS,prs_best_validation_EAS)

if(trait %in% c("Breast","Prostate")){
  confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
}else{
  confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
}

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]

roc_obj <- roc.binary(status = trait,
                      variable = "prs",
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

SL.result <- data.frame(method = "SL_Combined_Eur",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_NonEUR[!is.na(pheno_vad_NonEUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]

roc_obj <- roc.binary(status = trait,
                      variable = "prs",
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

SL.result <- data.frame(method = "SL_Combined_NonEur",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]

roc_obj <- roc.binary(status = trait,
                      variable = "prs",
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

SL.result <- data.frame(method = "SL_Combined_UNK",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]

roc_obj <- roc.binary(status = trait,
                      variable = "prs",
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

SL.result <- data.frame(method = "SL_Combined_SAS",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]

roc_obj <- roc.binary(status = trait,
                      variable = "prs",
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

SL.result <- data.frame(method = "SL_Combined_MIX",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]

roc_obj <- roc.binary(status = trait,
                      variable = "prs",
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

SL.result <- data.frame(method = "SL_Combined_AFR",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR.RData"))

## bootstrap the AUC to provide an approximate distribution 
if(trait %in% c("Prostate","CAD")){
  SL.result <- NA
}else{
  d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = "prs",
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  
  SL.result <- data.frame(method = "SL_Combined_EAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
}
save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS.RData"))

