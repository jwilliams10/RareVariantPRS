rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

STAARO_GeneCentric_Coding_Train_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/STAARO_GeneCentric_Coding_Train_PRS.csv")
STAARO_GeneCentric_Coding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/STAARO_GeneCentric_Coding_Tune_PRS.csv")
STAARO_GeneCentric_Coding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/STAARO_GeneCentric_Coding_Validation_PRS.csv")

STAARO_GeneCentric_Noncoding_Train_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/STAARO_GeneCentric_Noncoding_Train_PRS.csv")
STAARO_GeneCentric_Noncoding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/STAARO_GeneCentric_Noncoding_Tune_PRS.csv")
STAARO_GeneCentric_Noncoding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/STAARO_GeneCentric_Noncoding_Validation_PRS.csv")

STAARO_SlidingWindow_Train_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/STAARO_SlidingWindow_Train_PRS.csv")
STAARO_SlidingWindow_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/STAARO_SlidingWindow_Tune_PRS.csv")
STAARO_SlidingWindow_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/STAARO_SlidingWindow_Validation_PRS.csv")

colnames(STAARO_GeneCentric_Coding_Train_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Train_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Train_PRS)[2:ncol(STAARO_GeneCentric_Coding_Train_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Train_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Train_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Train_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Train_PRS)]))
colnames(STAARO_SlidingWindow_Train_PRS) <- c(colnames(STAARO_SlidingWindow_Train_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Train_PRS)[2:ncol(STAARO_SlidingWindow_Train_PRS)]))
STAARO_Combined_Train <- cbind(STAARO_GeneCentric_Coding_Train_PRS,STAARO_GeneCentric_Noncoding_Train_PRS[,-1],STAARO_SlidingWindow_Train_PRS[,-1])

colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Coding_Tune_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Tune_PRS)]))
colnames(STAARO_SlidingWindow_Tune_PRS) <- c(colnames(STAARO_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Tune_PRS)[2:ncol(STAARO_SlidingWindow_Tune_PRS)]))
STAARO_Combined_Tune <- cbind(STAARO_GeneCentric_Coding_Tune_PRS,STAARO_GeneCentric_Noncoding_Tune_PRS[,-1],STAARO_SlidingWindow_Tune_PRS[,-1])

colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Coding_Validation_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Validation_PRS)]))
colnames(STAARO_SlidingWindow_Validation_PRS) <- c(colnames(STAARO_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Validation_PRS)[2:ncol(STAARO_SlidingWindow_Validation_PRS)]))
STAARO_Combined_Validation <- cbind(STAARO_GeneCentric_Coding_Validation_PRS,STAARO_GeneCentric_Noncoding_Validation_PRS[,-1],STAARO_SlidingWindow_Validation_PRS[,-1])


Burden_GeneCentric_Coding_Train_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Burden_GeneCentric_Coding_Train_PRS.csv")
Burden_GeneCentric_Coding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Burden_GeneCentric_Coding_Tune_PRS.csv")
Burden_GeneCentric_Coding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Burden_GeneCentric_Coding_Validation_PRS.csv")

Burden_GeneCentric_Noncoding_Train_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Burden_GeneCentric_Noncoding_Train_PRS.csv")
Burden_GeneCentric_Noncoding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Burden_GeneCentric_Noncoding_Tune_PRS.csv")
Burden_GeneCentric_Noncoding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Burden_GeneCentric_Noncoding_Validation_PRS.csv")

Burden_SlidingWindow_Train_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Burden_SlidingWindow_Train_PRS.csv")
Burden_SlidingWindow_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Burden_SlidingWindow_Tune_PRS.csv")
Burden_SlidingWindow_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Burden_SlidingWindow_Validation_PRS.csv")

colnames(Burden_GeneCentric_Coding_Train_PRS) <- c(colnames(Burden_GeneCentric_Coding_Train_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Train_PRS)[2:ncol(Burden_GeneCentric_Coding_Train_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Train_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Train_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Train_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Train_PRS)]))
colnames(Burden_SlidingWindow_Train_PRS) <- c(colnames(Burden_SlidingWindow_Train_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Train_PRS)[2:ncol(Burden_SlidingWindow_Train_PRS)]))
Burden_Combined_Train <- cbind(Burden_GeneCentric_Coding_Train_PRS,Burden_GeneCentric_Noncoding_Train_PRS[,-1],Burden_SlidingWindow_Train_PRS[,-1])

colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Tune_PRS)[2:ncol(Burden_GeneCentric_Coding_Tune_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Tune_PRS)]))
colnames(Burden_SlidingWindow_Tune_PRS) <- c(colnames(Burden_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Tune_PRS)[2:ncol(Burden_SlidingWindow_Tune_PRS)]))
Burden_Combined_Tune <- cbind(Burden_GeneCentric_Coding_Tune_PRS,Burden_GeneCentric_Noncoding_Tune_PRS[,-1],Burden_SlidingWindow_Tune_PRS[,-1])

colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Validation_PRS)[2:ncol(Burden_GeneCentric_Coding_Validation_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Validation_PRS)]))
colnames(Burden_SlidingWindow_Validation_PRS) <- c(colnames(Burden_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Validation_PRS)[2:ncol(Burden_SlidingWindow_Validation_PRS)]))
Burden_Combined_Validation <- cbind(Burden_GeneCentric_Coding_Validation_PRS,Burden_GeneCentric_Noncoding_Validation_PRS[,-1],Burden_SlidingWindow_Validation_PRS[,-1])


# STAARO_Combined_Tune <- STAARO_Combined_Tune[,-c(16,31,46)]
# STAARO_Combined_Validation <- STAARO_Combined_Validation[,-c(16,31,46)]


## Drop Correlated Values
mtx <- cor(STAARO_Combined_Tune)
drop <- findCorrelation(mtx,cutoff=0.98)
drop <- names(STAARO_Combined_Tune)[drop]

STAARO_Combined_Train = STAARO_Combined_Train %>% 
  select(-all_of(drop))
STAARO_Combined_Tune = STAARO_Combined_Tune %>% 
  select(-all_of(drop))
STAARO_Combined_Validation = STAARO_Combined_Validation %>% 
  select(-all_of(drop))

mtx <- cor(Burden_Combined_Tune)
drop <- findCorrelation(mtx,cutoff=0.98)
drop <- names(Burden_Combined_Tune)[drop]

Burden_Combined_Train = Burden_Combined_Train %>% 
  select(-all_of(drop))
Burden_Combined_Tune = Burden_Combined_Tune %>% 
  select(-all_of(drop))
Burden_Combined_Validation = Burden_Combined_Validation %>% 
  select(-all_of(drop))

## Pull in Phenotypes/Covariates 
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt")
colnames(pheno_train) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

pheno_train_STAARO <- left_join(pheno_train,STAARO_Combined_Train,by = "IID")
STAARO_Combined_Train <- pheno_train_STAARO[,-c(1:16)]

pheno_train_Burden <- left_join(pheno_train,Burden_Combined_Train,by = "IID")
Burden_Combined_Train <- pheno_train_Burden[,-c(1:16)]

pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
colnames(pheno_tuning) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

pheno_tuning_STAARO <- left_join(pheno_tuning,STAARO_Combined_Tune,by = "IID")
STAARO_Combined_Tune <- pheno_tuning_STAARO[,-c(1:16)]

pheno_tuning_Burden <- left_join(pheno_tuning,Burden_Combined_Tune,by = "IID")
Burden_Combined_Tune <- pheno_tuning_Burden[,-c(1:16)]

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

pheno_vad_STAARO <- left_join(pheno_vad,STAARO_Combined_Validation,by = "IID")
STAARO_Combined_Validation <- pheno_vad_STAARO[,-c(1:16)]

pheno_vad_Burden <- left_join(pheno_vad,Burden_Combined_Validation,by = "IID")
Burden_Combined_Validation <- pheno_vad_Burden[,-c(1:16)]

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_train_STAARO)
y_train_STAARO <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_train_Burden)
y_train_Burden <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tuning_STAARO)
y_tune_STAARO <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tuning_Burden)
y_tune_Burden <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad_STAARO)
y_vad_STAARO <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad_Burden)
y_vad_Burden <- model.null$residual


arrayid <- as.numeric(commandArgs(TRUE)[1])

if(arrayid == 1){
  ##### SL 
  
  SL.library <- c(
    "SL.glmnet",
    "SL.ridge",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = y_tune_STAARO, X = STAARO_Combined_Tune, family = gaussian(),
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library)
  cvsl <- CV.SuperLearner(Y = y_tune_STAARO, X = STAARO_Combined_Tune, family = gaussian(),
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library)
  
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
  
  save(final_coefs,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/final_coefs_All_STAARO.RData")
  
  a <- predict(sl, STAARO_Combined_Validation, onlySL = FALSE)
  
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
  
  # prs_best_validation <- predict(sl, STAARO_Combined_Validation, onlySL = FALSE)[[1]]
  
  model <- lm(y_vad_STAARO~prs_best_validation)
  r2 <- summary(model)$r.square
  
  prs_best_validation <- data.frame(IID = pheno_vad_STAARO$IID,prs = prs_best_validation)
  
  
  a <- predict(sl, STAARO_Combined_Train, onlySL = FALSE)

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
  prs_best_train <- data.frame(IID = pheno_train_STAARO$IID,prs = prs_best_train)
  
  
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
  
  write.table(prs_best_train,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_STAARO_Train_All.txt",sep = "\t",row.names = FALSE)
  write.table(prs_best_tune,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_STAARO_Tune_All.txt",sep = "\t",row.names = FALSE)
  write.table(prs_best_validation,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_STAARO_Validation_All.txt",sep = "\t",row.names = FALSE)
  
  data <- data.frame(y = y_vad_STAARO, x = prs_best_validation$prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 100000)
  
  ci_result <- boot.ci(boot_r2, type = "bca")
  SL.result <- data.frame(method = "SL_STAARO",
                          r2 = r2,
                          r2_low = ci_result$bca[4],
                          r2_high = ci_result$bca[5]
  )
  
  save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/sl_result_All_STAARO.RData") 
}else{
  ##### SL 
  
  SL.library <- c(
    "SL.glmnet",
    "SL.ridge",
    "SL.glm",
    "SL.mean"
  )
  sl <- SuperLearner(Y = y_tune_Burden, X = Burden_Combined_Tune, family = gaussian(),
                     # For a real analysis we would use V = 10.
                     # V = 3,
                     SL.library = SL.library)
  cvsl <- CV.SuperLearner(Y = y_tune_Burden, X = Burden_Combined_Tune, family = gaussian(),
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library)
  
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
  
  save(final_coefs,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/final_coefs_All_Burden.RData")
  
  a <- predict(sl, Burden_Combined_Validation, onlySL = FALSE)
  
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
  
  # prs_best_validation <- predict(sl, Burden_Combined_Validation, onlySL = FALSE)[[1]]
  
  model <- lm(y_vad_Burden~prs_best_validation)
  r2 <- summary(model)$r.square
  
  prs_best_validation <- data.frame(IID = pheno_vad_Burden$IID,prs = prs_best_validation)
  
  
  a <- predict(sl, Burden_Combined_Train, onlySL = FALSE)
  
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
  prs_best_train <- data.frame(IID = pheno_train_Burden$IID,prs = prs_best_train)
  
  
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
  
  write.table(prs_best_train,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_Burden_Train_All.txt",sep = "\t",row.names = FALSE)
  write.table(prs_best_tune,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_Burden_Tune_All.txt",sep = "\t",row.names = FALSE)
  write.table(prs_best_validation,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/Best_All_Burden_Validation_All.txt",sep = "\t",row.names = FALSE)
  
  data <- data.frame(y = y_vad_Burden, x = prs_best_validation$prs)
  R2Boot <- function(data,indices){
    boot_data <- data[indices, ]
    model <- lm(y ~ x, data = boot_data)
    result <- summary(model)$r.square
    return(c(result))
  }
  boot_r2 <- boot(data = data, statistic = R2Boot, R = 100000)
  
  ci_result <- boot.ci(boot_r2, type = "bca")
  SL.result <- data.frame(method = "SL_Burden",
                          r2 = r2,
                          r2_low = ci_result$bca[4],
                          r2_high = ci_result$bca[5]
  )
  
  save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/sl_result_All_Burden.RData")
}
