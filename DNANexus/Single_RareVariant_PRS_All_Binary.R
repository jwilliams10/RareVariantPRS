rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(RISCA)
library(boot)
library(stringr)

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Single_RareVariant_PRS_All_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Single_RareVariant_PRS_All_Binary.sh -icmd="bash Single_RareVariant_PRS_All_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/ --priority high --instance-type mem1_ssd1_v2_x4

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  Noncoding_Burden <- read.csv(paste0(trait,"_Noncoding_Burden_PRS.csv"))
  Noncoding_STAARO <- read.csv(paste0(trait,"_Noncoding_STAARO_PRS.csv"))
  
  Coding_Burden <- read.csv(paste0(trait,"_Coding_Burden_PRS.csv"))
  Coding_STAARO <- read.csv(paste0(trait,"_Coding_STAARO_PRS.csv"))
  
  ncRNA_Burden <- read.csv(paste0(trait,"_ncRNA_Burden_PRS.csv"))
  ncRNA_STAARO <- read.csv(paste0(trait,"_ncRNA_STAARO_PRS.csv"))
  
  Tune_Null_Model <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
  Validation_Null_Model <- get(load(paste0(trait,"_Validation_Null_Model.RData")))
  
  file.remove(paste0(trait,"_Noncoding_Burden_PRS.csv"))
  file.remove(paste0(trait,"_Noncoding_STAARO_PRS.csv"))
  file.remove(paste0(trait,"_Coding_Burden_PRS.csv"))
  file.remove(paste0(trait,"_Coding_STAARO_PRS.csv"))
  file.remove(paste0(trait,"_ncRNA_Burden_PRS.csv"))
  file.remove(paste0(trait,"_ncRNA_STAARO_PRS.csv"))
  file.remove(paste0(trait,"_Tune_Null_Model.RData"))
  file.remove(paste0(trait,"_Validation_Null_Model.RData"))
  
  
  for(i in 2:16){
    Noncoding_Burden[,i] <- Noncoding_Burden[,i] + ncRNA_Burden[,i]
    Noncoding_STAARO[,i] <- Noncoding_STAARO[,i] + ncRNA_STAARO[,i]
  }
  rm(ncRNA_Burden);rm(ncRNA_STAARO)
  
  
  
  Noncoding_Burden_Tune <- Noncoding_Burden[Noncoding_Burden$IID %in% Tune_Null_Model$id_include,]
  Noncoding_Burden_Validation <- Noncoding_Burden[Noncoding_Burden$IID %in% Validation_Null_Model$id_include,]
  
  Noncoding_STAARO_Tune <- Noncoding_STAARO[Noncoding_STAARO$IID %in% Tune_Null_Model$id_include,]
  Noncoding_STAARO_Validation <- Noncoding_STAARO[Noncoding_STAARO$IID %in% Validation_Null_Model$id_include,]
  
  
  Coding_Burden_Tune <- Coding_Burden[Coding_Burden$IID %in% Tune_Null_Model$id_include,]
  Coding_Burden_Validation <- Coding_Burden[Coding_Burden$IID %in% Validation_Null_Model$id_include,]
  
  Coding_STAARO_Tune <- Coding_STAARO[Coding_STAARO$IID %in% Tune_Null_Model$id_include,]
  Coding_STAARO_Validation <- Coding_STAARO[Coding_STAARO$IID %in% Validation_Null_Model$id_include,]
  
  
  
  Burden_Combined_Tune <- cbind(Noncoding_Burden_Tune,Coding_Burden_Tune[,-1])
  colnames(Burden_Combined_Tune) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  STAARO_Combined_Tune <- cbind(Noncoding_STAARO_Tune,Coding_STAARO_Tune[,-1])
  colnames(STAARO_Combined_Tune) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  
  Burden_Combined_Validation <- cbind(Noncoding_Burden_Validation,Coding_Burden_Validation[,-1])
  colnames(Burden_Combined_Validation) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  STAARO_Combined_Validation <- cbind(Noncoding_STAARO_Validation,Coding_STAARO_Validation[,-1])
  colnames(STAARO_Combined_Validation) <- c("IID",paste0("Noncoding_Threshold",1:15),paste0("Coding_Threshold",1:15))
  
  
  
  
  drop <- caret::findLinearCombos(STAARO_Combined_Tune)$remove
  drop <- names(STAARO_Combined_Tune)[drop]
  
  STAARO_Combined_Tune = STAARO_Combined_Tune %>% 
    select(-all_of(drop))
  STAARO_Combined_Validation = STAARO_Combined_Validation %>% 
    select(-all_of(drop))
  
  ## Drop Correlated Values
  mtx <- cor(STAARO_Combined_Tune)
  drop <- findCorrelation(mtx,cutoff=0.98)
  drop <- names(STAARO_Combined_Tune)[drop]
  
  STAARO_Combined_Tune = STAARO_Combined_Tune %>% 
    select(-all_of(drop))
  STAARO_Combined_Validation = STAARO_Combined_Validation %>% 
    select(-all_of(drop))
  
  drop <- caret::findLinearCombos(Burden_Combined_Tune)$remove
  drop <- names(Burden_Combined_Tune)[drop]
  
  Burden_Combined_Tune = Burden_Combined_Tune %>% 
    select(-all_of(drop))
  Burden_Combined_Validation = Burden_Combined_Validation %>% 
    select(-all_of(drop))
  
  mtx <- cor(Burden_Combined_Tune)
  drop <- findCorrelation(mtx,cutoff=0.98)
  drop <- names(Burden_Combined_Tune)[drop]
  
  Burden_Combined_Tune = Burden_Combined_Tune %>% 
    select(-all_of(drop))
  Burden_Combined_Validation = Burden_Combined_Validation %>% 
    select(-all_of(drop))
  
  
  
  
  
  ## Pull in Phenotypes/Covariates 
  pheno_tuning <- read.delim("All_Tune.txt")
  
  pheno_tuning_STAARO <- left_join(pheno_tuning,STAARO_Combined_Tune,by = "IID")
  pheno_tuning_STAARO <- pheno_tuning_STAARO[!is.na(pheno_tuning_STAARO[,trait]),]
  STAARO_Combined_Tune <- pheno_tuning_STAARO[,-c(1:26)]
  
  pheno_tuning_Burden <- left_join(pheno_tuning,Burden_Combined_Tune,by = "IID")
  pheno_tuning_Burden <- pheno_tuning_Burden[!is.na(pheno_tuning_Burden[,trait]),]
  Burden_Combined_Tune <- pheno_tuning_Burden[,-c(1:26)]
  
  pheno_vad <- read.delim("All_Validation.txt")
  
  pheno_vad_STAARO <- left_join(pheno_vad,STAARO_Combined_Validation,by = "IID")
  pheno_vad_STAARO <- pheno_vad_STAARO[!is.na(pheno_vad_STAARO[,trait]),]
  STAARO_Combined_Validation <- pheno_vad_STAARO[,-c(1:26)]
  
  pheno_vad_Burden <- left_join(pheno_vad,Burden_Combined_Validation,by = "IID")
  pheno_vad_Burden <- pheno_vad_Burden[!is.na(pheno_vad_Burden[,trait]),]
  Burden_Combined_Validation <- pheno_vad_Burden[,-c(1:26)]
  
  
  
  if(trait %in% c("Breast","Prostate")){
    confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
  }else{
    confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
  }
  
  
  
  
  
  
  SL.library <- c(
    "SL.glmnet",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = pheno_tuning_STAARO[,trait], X = STAARO_Combined_Tune, family = binomial(), method = "method.AUC",
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library,cvControl = list(V = 20,stratifyCV = TRUE))
  cvsl <- CV.SuperLearner(Y = pheno_tuning_STAARO[,trait], X = STAARO_Combined_Tune, family = binomial(), method = "method.AUC",
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library,cvControl = list(V = 20,stratifyCV = TRUE))
  
  best_algorithm <- summary(cvsl)$Table$Algorithm[which.max(summary(cvsl)$Table$Ave)]
  
  ### Extract Coefficients
  #algorithm weight
  alg_weights <- sl$coef
  #glmnet
  glmnet_obj <- sl$fitLibrary$SL.glmnet$object
  best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
  glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
  #glm 
  glm_coefs <- sl$fitLibrary$SL.glm_All$object$coefficients
  #mean 
  mean_coefs <- sl$fitLibrary$SL.mean_All$object
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    final_coefs <- glmnet_coefs
    final_coefs = final_coefs[2:nrow(final_coefs),]
  }else if(best_algorithm == "SL.glm_All"){
    #final
    final_coefs <- glm_coefs
    final_coefs = final_coefs[2:length(final_coefs)]
  }else{
    #final
    final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * glm_coefs + alg_weights[3] * mean_coefs 
    final_coefs = final_coefs[2:nrow(final_coefs),]
  }

  #remove weight 0 coefficients
  final_coefs = final_coefs[final_coefs!=0]
  
  save(final_coefs,file = paste0(trait,"_final_coefs_All_STAARO.RData"))
  
  if(best_algorithm == "SL.glmnet_All"){
    a <- predict(sl, STAARO_Combined_Validation,onlySL = FALSE)
    prs_best_validation_glmnet <- a$library.predict[,1]
    #final
    prs_best_validation <- prs_best_validation_glmnet
  }else if(best_algorithm == "SL.glm_All"){
    a <- predict(sl, STAARO_Combined_Validation,onlySL = FALSE)
    prs_best_validation_glm <- a$library.predict[,3]
    #final
    prs_best_validation <- prs_best_validation_glm
  }else{
    a <- predict(sl, STAARO_Combined_Validation,onlySL = TRUE)
    prs_best_validation_sl <- a$pred
    #final
    prs_best_validation <- prs_best_validation_sl
  }
  
  prs_best_validation <- data.frame(IID = pheno_vad_STAARO$IID,prs = prs_best_validation)
  
  a <- predict(sl, STAARO_Combined_Tune, onlySL = FALSE)
  
  prs_best_tune_sl <- a$pred
  prs_best_tune_glmnet <- a$library.predict[,1]
  prs_best_tune_glm <- a$library.predict[,2]
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    prs_best_tune <- prs_best_tune_glmnet
  }else if(best_algorithm == "SL.glm_All"){
    #final
    prs_best_tune <- prs_best_tune_glm
  }else{
    #final
    prs_best_tune <- prs_best_tune_sl
  }
  prs_best_tune <- data.frame(IID = pheno_tuning_STAARO$IID,prs = prs_best_tune)
  
  pheno_vad_STAARO <- left_join(pheno_vad_STAARO,prs_best_validation)
  pheno_tuning_STAARO <- left_join(pheno_tuning_STAARO,prs_best_tune)
  
  prs_columns <- c(which(str_detect(colnames(pheno_tuning_STAARO),"Coding_Threshold")),which(str_detect(colnames(pheno_tuning_STAARO),"Noncoding_Threshold")),which(str_detect(colnames(pheno_tuning_STAARO),"prs")))
  
  auc_tune <- vector()
  for(i in 1:length(prs_columns)){
    
    auc_tune[i] <- roc.binary(status = trait,
                              variable = colnames(pheno_tuning_STAARO)[prs_columns[i]],
                              confounders = confounders,
                              data = pheno_tuning_STAARO[!is.na(pheno_tuning_STAARO[,trait]),],
                              precision=seq(0.05,0.95, by=0.05))$auc
  }
  prs_best_tune <- data.frame(IID = pheno_tuning_STAARO$IID,prs = pheno_tuning_STAARO[,colnames(pheno_tuning_STAARO)[prs_columns[which.max(auc_tune)]]])
  pheno_tuning_STAARO <- pheno_tuning_STAARO[,colnames(pheno_tuning_STAARO) != "prs"]
  pheno_tuning_STAARO <- left_join(pheno_tuning_STAARO,prs_best_tune)
  
  prs_best_validation <- data.frame(IID = pheno_vad_STAARO$IID,prs = pheno_vad_STAARO[,colnames(pheno_tuning_STAARO)[prs_columns[which.max(auc_tune)]]])
  pheno_vad_STAARO <- pheno_vad_STAARO[,colnames(pheno_vad_STAARO) != "prs"]
  pheno_vad_STAARO <- left_join(pheno_vad_STAARO,prs_best_validation)
  
  write.table(prs_best_tune,file=paste0(trait,"_Best_All_STAARO_Tune_All.txt"),sep = "\t",row.names = FALSE)
  write.table(prs_best_validation,file=paste0(trait,"_Best_All_STAARO_Validation_All.txt"),sep = "\t",row.names = FALSE)
  
  load("all_phenotypes.RData")
  
  pheno_vad_EUR <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_vad_NonEur <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_vad_UNK <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_vad_SAS <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_vad_MIX <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_vad_AFR <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_vad_EAS <- pheno_vad_STAARO[pheno_vad_STAARO$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  ## bootstrap the R2 to provide an approximate distribution 
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
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  SL.result <- data.frame(method = "SL_Combined_EUR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_Eur_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_NonEUR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_NonEur_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  SL.result <- data.frame(method = "SL_Combined_UNK",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_UNK_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_SAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_SAS_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_MIX",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_MIX_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_AFR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_AFR_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_EAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_EAS_STAARO.RData"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  SL.library <- c(
    "SL.glmnet",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = pheno_tuning_Burden[,trait], X = Burden_Combined_Tune, family = binomial(), method = "method.AUC",
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library,cvControl = list(V = 20,stratifyCV = TRUE))
  cvsl <- CV.SuperLearner(Y = pheno_tuning_Burden[,trait], X = Burden_Combined_Tune, family = binomial(), method = "method.AUC",
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library,cvControl = list(V = 20,stratifyCV = TRUE))
  
  best_algorithm <- summary(cvsl)$Table$Algorithm[which.max(summary(cvsl)$Table$Ave)]
  
  ### Extract Coefficients
  #algorithm weight
  alg_weights <- sl$coef
  #glmnet
  glmnet_obj <- sl$fitLibrary$SL.glmnet$object
  best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
  glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
  #glm 
  glm_coefs <- sl$fitLibrary$SL.glm_All$object$coefficients
  #mean 
  mean_coefs <- sl$fitLibrary$SL.mean_All$object
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    final_coefs <- glmnet_coefs
    final_coefs = final_coefs[2:nrow(final_coefs),]
  }else if(best_algorithm == "SL.glm_All"){
    #final
    final_coefs <- glm_coefs
    final_coefs = final_coefs[2:length(final_coefs)]
  }else{
    #final
    final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * glm_coefs + alg_weights[3] * mean_coefs 
    final_coefs = final_coefs[2:nrow(final_coefs),]
  }
  
  #remove weight 0 coefficients
  final_coefs = final_coefs[final_coefs!=0]
  
  save(final_coefs,file = paste0(trait,"_final_coefs_All_Burden.RData"))
  
  if(best_algorithm == "SL.glmnet_All"){
    a <- predict(sl, Burden_Combined_Validation,onlySL = FALSE)
    prs_best_validation_glmnet <- a$library.predict[,1]
    #final
    prs_best_validation <- prs_best_validation_glmnet
  }else if(best_algorithm == "SL.glm_All"){
    a <- predict(sl, Burden_Combined_Validation,onlySL = FALSE)
    prs_best_validation_glm <- a$library.predict[,3]
    #final
    prs_best_validation <- prs_best_validation_glm
  }else{
    a <- predict(sl, Burden_Combined_Validation,onlySL = TRUE)
    prs_best_validation_sl <- a$pred
    #final
    prs_best_validation <- prs_best_validation_sl
  }
  
  prs_best_validation <- data.frame(IID = pheno_vad_Burden$IID,prs = prs_best_validation)
  
  a <- predict(sl, Burden_Combined_Tune, onlySL = FALSE)
  
  prs_best_tune_sl <- a$pred
  prs_best_tune_glmnet <- a$library.predict[,1]
  prs_best_tune_glm <- a$library.predict[,2]
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    prs_best_tune <- prs_best_tune_glmnet
  }else if(best_algorithm == "SL.glm_All"){
    #final
    prs_best_tune <- prs_best_tune_glm
  }else{
    #final
    prs_best_tune <- prs_best_tune_sl
  }
  prs_best_tune <- data.frame(IID = pheno_tuning_Burden$IID,prs = prs_best_tune)
  
  pheno_vad_Burden <- left_join(pheno_vad_Burden,prs_best_validation)
  pheno_tuning_Burden <- left_join(pheno_tuning_Burden,prs_best_tune)
  
  prs_columns <- c(which(str_detect(colnames(pheno_tuning_Burden),"Coding_Threshold")),which(str_detect(colnames(pheno_tuning_Burden),"Noncoding_Threshold")),which(str_detect(colnames(pheno_tuning_Burden),"prs")))
  
  auc_tune <- vector()
  for(i in 1:length(prs_columns)){
    
    auc_tune[i] <- roc.binary(status = trait,
                              variable = colnames(pheno_tuning_Burden)[prs_columns[i]],
                              confounders = confounders,
                              data = pheno_tuning_Burden[!is.na(pheno_tuning_Burden[,trait]),],
                              precision=seq(0.05,0.95, by=0.05))$auc
  }
  prs_best_tune <- data.frame(IID = pheno_tuning_Burden$IID,prs = pheno_tuning_Burden[,colnames(pheno_tuning_Burden)[prs_columns[which.max(auc_tune)]]])
  pheno_tuning_Burden <- pheno_tuning_Burden[,colnames(pheno_tuning_Burden) != "prs"]
  pheno_tuning_Burden <- left_join(pheno_tuning_Burden,prs_best_tune)
  
  prs_best_validation <- data.frame(IID = pheno_vad_Burden$IID,prs = pheno_vad_Burden[,colnames(pheno_tuning_Burden)[prs_columns[which.max(auc_tune)]]])
  pheno_vad_Burden <- pheno_vad_Burden[,colnames(pheno_vad_Burden) != "prs"]
  pheno_vad_Burden <- left_join(pheno_vad_Burden,prs_best_validation)
  
  write.table(prs_best_tune,file=paste0(trait,"_Best_All_Burden_Tune_All.txt"),sep = "\t",row.names = FALSE)
  write.table(prs_best_validation,file=paste0(trait,"_Best_All_Burden_Validation_All.txt"),sep = "\t",row.names = FALSE)
  
  load("all_phenotypes.RData")
  
  pheno_vad_EUR <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_vad_NonEur <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_vad_UNK <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_vad_SAS <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_vad_MIX <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_vad_AFR <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_vad_EAS <- pheno_vad_Burden[pheno_vad_Burden$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_EUR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_Eur_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_NonEUR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_NonEur_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_UNK",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_UNK_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_SAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_SAS_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_MIX",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_MIX_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_AFR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_AFR_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","prs")]
  roc_obj <- roc.binary(status = trait,
                        variable = "prs",
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  if(AUC == 0.5){
    ci_result <- data.frame(percent = c(0,0,0,.5,.5))
  }else{
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
  }
  SL.result <- data.frame(method = "SL_Combined_EAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  save(SL.result,file = paste0(trait,"_sl_result_All_EAS_Burden.RData"))
  
  
}

system("rm All_Tune.txt")
system("rm All_Validation.txt")
system("rm all_phenotypes.RData")