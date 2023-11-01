rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))
obj_nullmodel_sub <- get(load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Tune_Null_Model1.RData"))

ids_tune <- obj_nullmodel_full$id_include[obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include]

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))
obj_nullmodel_sub <- get(load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Validation_Null_Model1.RData"))

ids_validation <- obj_nullmodel_full$id_include[obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include]


load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/tune_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
STAARO_GeneCentric_Coding_Tune_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Coding_Tune_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/validation_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
STAARO_GeneCentric_Coding_Validation_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Coding_Validation_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/tune_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
STAARO_GeneCentric_Noncoding_Tune_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Noncoding_Tune_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/validation_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
STAARO_GeneCentric_Noncoding_Validation_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Noncoding_Validation_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/tune_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
STAARO_SlidingWindow_Tune_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_SlidingWindow_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_SlidingWindow_Tune_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_SlidingWindow_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/validation_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
STAARO_SlidingWindow_Validation_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_SlidingWindow_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_SlidingWindow_Validation_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_SlidingWindow_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Coding_Tune_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Tune_PRS)]))
colnames(STAARO_SlidingWindow_Tune_PRS) <- c(colnames(STAARO_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Tune_PRS)[2:ncol(STAARO_SlidingWindow_Tune_PRS)]))
STAARO_Combined_Tune <- cbind(STAARO_GeneCentric_Coding_Tune_PRS,STAARO_GeneCentric_Noncoding_Tune_PRS[,-1],STAARO_SlidingWindow_Tune_PRS[,-1])

colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Coding_Validation_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Validation_PRS)]))
colnames(STAARO_SlidingWindow_Validation_PRS) <- c(colnames(STAARO_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Validation_PRS)[2:ncol(STAARO_SlidingWindow_Validation_PRS)]))
STAARO_Combined_Validation <- cbind(STAARO_GeneCentric_Coding_Validation_PRS,STAARO_GeneCentric_Noncoding_Validation_PRS[,-1],STAARO_SlidingWindow_Validation_PRS[,-1])

colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Tune_PRS)[2:ncol(Burden_GeneCentric_Coding_Tune_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Tune_PRS)]))
colnames(Burden_SlidingWindow_Tune_PRS) <- c(colnames(Burden_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Tune_PRS)[2:ncol(Burden_SlidingWindow_Tune_PRS)]))
Burden_Combined_Tune <- cbind(Burden_GeneCentric_Coding_Tune_PRS,Burden_GeneCentric_Noncoding_Tune_PRS[,-1],Burden_SlidingWindow_Tune_PRS[,-1])

colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Validation_PRS)[2:ncol(Burden_GeneCentric_Coding_Validation_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Validation_PRS)]))
colnames(Burden_SlidingWindow_Validation_PRS) <- c(colnames(Burden_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Validation_PRS)[2:ncol(Burden_SlidingWindow_Validation_PRS)]))
Burden_Combined_Validation <- cbind(Burden_GeneCentric_Coding_Validation_PRS,Burden_GeneCentric_Noncoding_Validation_PRS[,-1],Burden_SlidingWindow_Validation_PRS[,-1])


## Drop Correlated Values
drop <- caret::findLinearCombos(STAARO_Combined_Tune)$remove
drop <- names(STAARO_Combined_Tune)[drop]

STAARO_Combined_Tune = STAARO_Combined_Tune %>% 
  select(-all_of(drop))
STAARO_Combined_Validation = STAARO_Combined_Validation %>% 
  select(-all_of(drop))

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
pheno_tuning <- Y_tune[[i]]
colnames(pheno_tuning) <- c("IID","Y")

pheno_tuning_STAARO <- left_join(pheno_tuning,STAARO_Combined_Tune,by = "IID")
STAARO_Combined_Tune <- pheno_tuning_STAARO[,-c(1:2)]

pheno_tuning_Burden <- left_join(pheno_tuning,Burden_Combined_Tune,by = "IID")
Burden_Combined_Tune <- pheno_tuning_Burden[,-c(1:2)]

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")

pheno_vad_STAARO <- left_join(pheno_vad,STAARO_Combined_Validation,by = "IID")
STAARO_Combined_Validation <- pheno_vad_STAARO[,-c(1:2)]

pheno_vad_Burden <- left_join(pheno_vad,Burden_Combined_Validation,by = "IID")
Burden_Combined_Validation <- pheno_vad_Burden[,-c(1:2)]

model.null <- lm(Y~1,data=pheno_tuning_STAARO)
y_tune_STAARO <- model.null$residual

model.null <- lm(Y~1,data=pheno_tuning_Burden)
y_tune_Burden <- model.null$residual

model.null <- lm(Y~1,data=pheno_vad_STAARO)
y_vad_STAARO <- model.null$residual

model.null <- lm(Y~1,data=pheno_vad_Burden)
y_vad_Burden <- model.null$residual

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

model <- lm(y_vad_STAARO~prs_best_validation)
r2 <- summary(model)$r.square

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

write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_STAARO_Tune_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_STAARO_Validation_All",i,".txt"),sep = "\t",row.names = FALSE)

data <- data.frame(y = y_vad_STAARO, x = prs_best_validation$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO",i,".RData"))

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

model <- lm(y_vad_Burden~prs_best_validation)
r2 <- summary(model)$r.square

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

write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_Burden_Tune_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_Burden_Validation_All",i,".txt"),sep = "\t",row.names = FALSE)

data <- data.frame(y = y_vad_Burden, x = prs_best_validation$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden",i,".RData"))
  
