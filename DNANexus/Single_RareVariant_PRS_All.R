rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(stringr)

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Single_RareVariant_PRS_All.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Single_RareVariant_PRS_All.sh -icmd="bash Single_RareVariant_PRS_All.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/ --priority low --instance-type mem1_ssd1_v2_x4

trait <- "BMI"

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
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
  
  
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_STAARO)
  y_tune_STAARO <- model.null$residual
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning_Burden)
  y_tune_Burden <- model.null$residual
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_STAARO)
  y_vad_STAARO <- model.null$residual
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_Burden)
  y_vad_Burden <- model.null$residual
  
  
  pheno_tuning_STAARO$y_tune <- NA
  pheno_tuning_STAARO$y_tune[!is.na(pheno_tuning_STAARO[,trait])] <- y_tune_STAARO
  
  pheno_tuning_Burden$y_tune <- NA
  pheno_tuning_Burden$y_tune[!is.na(pheno_tuning_Burden[,trait])] <- y_tune_STAARO
  
  
  
  
  
  
  
  
  
  SL.library <- c(
    "SL.glmnet",
    "SL.ridge",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = y_tune_STAARO, X = STAARO_Combined_Tune, family = gaussian(),
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library,cvControl = list(V = 20))
  cvsl <- CV.SuperLearner(Y = y_tune_STAARO, X = STAARO_Combined_Tune, family = gaussian(),
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library,cvControl = list(V = 20))
  
  best_algorithm <- summary(cvsl)$Table$Algorithm[which.min(summary(cvsl)$Table$Ave)]
  
  ### Extract Coefficients
  #algorithm weight
  alg_weights <- sl$coef
  #glmnet
  glmnet_obj <- sl$fitLibrary$SL.glmnet$object
  best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
  glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
  #ridge
  ridge_coefs <- sl$fitLibrary$SL.ridge_All$bestCoef
  #glm 
  glm_coefs <- sl$fitLibrary$SL.glm_All$object$coefficients
  #mean 
  mean_coefs <- sl$fitLibrary$SL.mean_All$object
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    final_coefs <- glmnet_coefs
  }else if(best_algorithm == "SL.ridge_All"){
    #final
    final_coefs <- ridge_coefs
  }else if(best_algorithm == "SL.glm_All"){
    #final
    final_coefs <- glm_coefs
  }else{
    #final
    final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * ridge_coefs + alg_weights[3] * glm_coefs + alg_weights[4] * mean_coefs 
  }
  
  #remove the intercept
  final_coefs = final_coefs[2:nrow(final_coefs),]
  #remove weight 0 coefficients
  final_coefs = final_coefs[final_coefs!=0]
  
  save(final_coefs,file = paste0(trait,"_final_coefs_All_STAARO.RData"))
  
  if(best_algorithm == "SL.glmnet_All"){
    a <- predict(sl, STAARO_Combined_Validation,onlySL = FALSE)
    prs_best_validation_glmnet <- a$library.predict[,1]
    #final
    prs_best_validation <- prs_best_validation_glmnet
  }else if(best_algorithm == "SL.ridge_All"){
    a <- predict(sl, STAARO_Combined_Validation,onlySL = FALSE)
    prs_best_validation_ridge <- a$library.predict[,2]
    #final
    prs_best_validation <- prs_best_validation_ridge
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
  prs_best_tune <- data.frame(IID = pheno_tuning_STAARO$IID,prs = prs_best_tune)
  
  pheno_vad_STAARO <- left_join(pheno_vad_STAARO,prs_best_validation)
  pheno_tuning_STAARO <- left_join(pheno_tuning_STAARO,prs_best_tune)
  
  prs_columns <- c(which(str_detect(colnames(pheno_tuning_STAARO),"Coding_Threshold")),which(str_detect(colnames(pheno_tuning_STAARO),"Noncoding_Threshold")),which(str_detect(colnames(pheno_tuning_STAARO),"prs")))
  
  r2_tune <- vector()
  for(i in 1:length(prs_columns)){
    r2_tune[i] <- summary(lm(as.formula(paste0("y_tune ~",colnames(pheno_tuning_STAARO)[prs_columns[i]])),data = pheno_tuning_STAARO))$r.squared
  }
  prs_best_tune <- data.frame(IID = pheno_tuning_STAARO$IID,prs = pheno_tuning_STAARO[,colnames(pheno_tuning_STAARO)[prs_columns[which.max(r2_tune)]]])
  
  prs_best_validation <- data.frame(IID = pheno_vad_STAARO$IID,prs = pheno_vad_STAARO[,colnames(pheno_tuning_STAARO)[prs_columns[which.max(r2_tune)]]])
  
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
  
  prs_best_validation_EUR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  prs_best_validation_NonEur <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  prs_best_validation_UNK <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  prs_best_validation_SAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  prs_best_validation_MIX <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  prs_best_validation_AFR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  prs_best_validation_EAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- prs_best_validation_EUR[!is.na(pheno_vad_EUR[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_Eur_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- prs_best_validation_NonEur[!is.na(pheno_vad_NonEur[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_NonEur_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- prs_best_validation_UNK[!is.na(pheno_vad_UNK[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_UNK_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- prs_best_validation_SAS[!is.na(pheno_vad_SAS[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_SAS_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- prs_best_validation_MIX[!is.na(pheno_vad_MIX[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_MIX_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- prs_best_validation_AFR[!is.na(pheno_vad_AFR[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_AFR_STAARO.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- prs_best_validation_EAS[!is.na(pheno_vad_EAS[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_EAS_STAARO.RData"))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  SL.library <- c(
    "SL.glmnet",
    "SL.ridge",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = y_tune_STAARO, X = Burden_Combined_Tune, family = gaussian(),
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library,cvControl = list(V = 20))
  cvsl <- CV.SuperLearner(Y = y_tune_Burden, X = Burden_Combined_Tune, family = gaussian(),
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library,cvControl = list(V = 20))
  
  best_algorithm <- summary(cvsl)$Table$Algorithm[which.min(summary(cvsl)$Table$Ave)]
  
  ### Extract Coefficients
  #algorithm weight
  alg_weights <- sl$coef
  #glmnet
  glmnet_obj <- sl$fitLibrary$SL.glmnet$object
  best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
  glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
  #ridge
  ridge_coefs <- sl$fitLibrary$SL.ridge_All$bestCoef
  #glm 
  glm_coefs <- sl$fitLibrary$SL.glm_All$object$coefficients
  #mean 
  mean_coefs <- sl$fitLibrary$SL.mean_All$object
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    final_coefs <- glmnet_coefs
  }else if(best_algorithm == "SL.ridge_All"){
    #final
    final_coefs <- ridge_coefs
  }else if(best_algorithm == "SL.glm_All"){
    #final
    final_coefs <- glm_coefs
  }else{
    #final
    final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * ridge_coefs + alg_weights[3] * glm_coefs + alg_weights[4] * mean_coefs 
  }
  
  #remove the intercept
  final_coefs = final_coefs[2:nrow(final_coefs),]
  #remove weight 0 coefficients
  final_coefs = final_coefs[final_coefs!=0]
  
  save(final_coefs,file = paste0(trait,"_final_coefs_All_Burden.RData"))
  
  if(best_algorithm == "SL.glmnet_All"){
    a <- predict(sl, Burden_Combined_Validation,onlySL = FALSE)
    prs_best_validation_glmnet <- a$library.predict[,1]
    #final
    prs_best_validation <- prs_best_validation_glmnet
  }else if(best_algorithm == "SL.ridge_All"){
    a <- predict(sl, Burden_Combined_Validation,onlySL = FALSE)
    prs_best_validation_ridge <- a$library.predict[,2]
    #final
    prs_best_validation <- prs_best_validation_ridge
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
  prs_best_tune <- data.frame(IID = pheno_tuning_Burden$IID,prs = prs_best_tune)
  
  pheno_vad_Burden <- left_join(pheno_vad_Burden,prs_best_validation)
  pheno_tuning_Burden <- left_join(pheno_tuning_Burden,prs_best_tune)
  
  prs_columns <- c(which(str_detect(colnames(pheno_tuning_Burden),"Coding_Threshold")),which(str_detect(colnames(pheno_tuning_Burden),"Noncoding_Threshold")),which(str_detect(colnames(pheno_tuning_Burden),"prs")))
  
  r2_tune <- vector()
  for(i in 1:length(prs_columns)){
    r2_tune[i] <- summary(lm(as.formula(paste0("y_tune ~",colnames(pheno_tuning_Burden)[prs_columns[i]])),data = pheno_tuning_Burden))$r.squared
  }
  prs_best_tune <- data.frame(IID = pheno_tuning_Burden$IID,prs = pheno_tuning_Burden[,colnames(pheno_tuning_Burden)[prs_columns[which.max(r2_tune)]]])
  
  prs_best_validation <- data.frame(IID = pheno_vad_Burden$IID,prs = pheno_vad_Burden[,colnames(pheno_tuning_Burden)[prs_columns[which.max(r2_tune)]]])
  
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
  
  prs_best_validation_EUR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  prs_best_validation_NonEur <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  prs_best_validation_UNK <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  prs_best_validation_SAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  prs_best_validation_MIX <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  prs_best_validation_AFR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  prs_best_validation_EAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EUR)
  prs <- prs_best_validation_EUR[!is.na(pheno_vad_EUR[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_Eur_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_NonEur)
  prs <- prs_best_validation_NonEur[!is.na(pheno_vad_NonEur[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_NonEur_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_UNK)
  prs <- prs_best_validation_UNK[!is.na(pheno_vad_UNK[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_UNK_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_SAS)
  prs <- prs_best_validation_SAS[!is.na(pheno_vad_SAS[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_SAS_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_MIX)
  prs <- prs_best_validation_MIX[!is.na(pheno_vad_MIX[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_MIX_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_AFR)
  prs <- prs_best_validation_AFR[!is.na(pheno_vad_AFR[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_AFR_Burden.RData"))
  
  ## bootstrap the R2 to provide an approximate distribution 
  model.vad.null  <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad_EAS)
  prs <- prs_best_validation_EAS[!is.na(pheno_vad_EAS[,trait]),2]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2 <- summary(model.vad.prs)$r.square
  
  data <- data.frame(y = model.vad.null$residual, x = prs)
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
  
  save(SL.result,file = paste0(trait,"_sl_result_All_EAS_Burden.RData"))
  
  
}

system("rm All_Tune.txt")
system("rm All_Validation.txt")
system("rm all_phenotypes.RData")