rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(stringr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))
obj_nullmodel_sub <- get(load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Tune_Null_Model1.RData"))

ids_tune <- obj_nullmodel_full$id_include[obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include]

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))
obj_nullmodel_sub <- get(load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Validation_Null_Model1.RData"))

ids_validation <- obj_nullmodel_full$id_include[obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include]


load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/tune_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
STAARO_GeneCentric_Coding_Tune_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Coding_Tune_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/validation_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
STAARO_GeneCentric_Coding_Validation_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Coding_Validation_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/tune_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
STAARO_GeneCentric_Noncoding_Tune_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Noncoding_Tune_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/validation_prs_mat",i,".RData"))
prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
STAARO_GeneCentric_Noncoding_Validation_PRS <- prs_mat[,c(1,2:16)]
colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
Burden_GeneCentric_Noncoding_Validation_PRS <- prs_mat[,c(1,17:31)]
colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

# load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/tune_prs_mat",i,".RData"))
# prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
# STAARO_SlidingWindow_Tune_PRS <- prs_mat[,c(1,2:16)]
# colnames(STAARO_SlidingWindow_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
# Burden_SlidingWindow_Tune_PRS <- prs_mat[,c(1,17:31)]
# colnames(Burden_SlidingWindow_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
# 
# load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/validation_prs_mat",i,".RData"))
# prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
# STAARO_SlidingWindow_Validation_PRS <- prs_mat[,c(1,2:16)]
# colnames(STAARO_SlidingWindow_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
# Burden_SlidingWindow_Validation_PRS <- prs_mat[,c(1,17:31)]
# colnames(Burden_SlidingWindow_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Coding_Tune_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Tune_PRS)]))
# colnames(STAARO_SlidingWindow_Tune_PRS) <- c(colnames(STAARO_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Tune_PRS)[2:ncol(STAARO_SlidingWindow_Tune_PRS)]))
STAARO_Combined_Tune <- cbind(STAARO_GeneCentric_Coding_Tune_PRS,STAARO_GeneCentric_Noncoding_Tune_PRS[,-1])#,STAARO_SlidingWindow_Tune_PRS[,-1])

colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(STAARO_GeneCentric_Coding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Coding_Validation_PRS)]))
colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[2:ncol(STAARO_GeneCentric_Noncoding_Validation_PRS)]))
# colnames(STAARO_SlidingWindow_Validation_PRS) <- c(colnames(STAARO_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(STAARO_SlidingWindow_Validation_PRS)[2:ncol(STAARO_SlidingWindow_Validation_PRS)]))
STAARO_Combined_Validation <- cbind(STAARO_GeneCentric_Coding_Validation_PRS,STAARO_GeneCentric_Noncoding_Validation_PRS[,-1])#,STAARO_SlidingWindow_Validation_PRS[,-1])

colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Coding_Tune_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Tune_PRS)[2:ncol(Burden_GeneCentric_Coding_Tune_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Tune_PRS)]))
# colnames(Burden_SlidingWindow_Tune_PRS) <- c(colnames(Burden_SlidingWindow_Tune_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Tune_PRS)[2:ncol(Burden_SlidingWindow_Tune_PRS)]))
Burden_Combined_Tune <- cbind(Burden_GeneCentric_Coding_Tune_PRS,Burden_GeneCentric_Noncoding_Tune_PRS[,-1])#,Burden_SlidingWindow_Tune_PRS[,-1])

colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Coding_Validation_PRS)[1],paste0("GeneCentric_Coding_",colnames(Burden_GeneCentric_Coding_Validation_PRS)[2:ncol(Burden_GeneCentric_Coding_Validation_PRS)]))
colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[1],paste0("GeneCentric_Noncoding_",colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[2:ncol(Burden_GeneCentric_Noncoding_Validation_PRS)]))
# colnames(Burden_SlidingWindow_Validation_PRS) <- c(colnames(Burden_SlidingWindow_Validation_PRS)[1],paste0("SlidingWindow_",colnames(Burden_SlidingWindow_Validation_PRS)[2:ncol(Burden_SlidingWindow_Validation_PRS)]))
Burden_Combined_Validation <- cbind(Burden_GeneCentric_Coding_Validation_PRS,Burden_GeneCentric_Noncoding_Validation_PRS[,-1])#,Burden_SlidingWindow_Validation_PRS[,-1])


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

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")

pheno_vad_STAARO <- left_join(pheno_vad,STAARO_Combined_Validation,by = "IID")
STAARO_Combined_Validation <- pheno_vad_STAARO[,-c(1:2)]

pheno_vad_Burden <- left_join(pheno_vad,Burden_Combined_Validation,by = "IID")
Burden_Combined_Validation <- pheno_vad_Burden[,-c(1:2)]

model.null <- lm(Y~1,data=pheno_tuning_STAARO)
y_tune_STAARO <- model.null$residual

pheno_tuning_STAARO$y_tune <- y_tune_STAARO

model.null <- lm(Y~1,data=pheno_tuning_Burden)
y_tune_Burden <- model.null$residual

pheno_tuning_Burden$y_tune <- y_tune_Burden

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
                   SL.library = SL.library,cvControl = list(V = 20))
cvsl <- CV.SuperLearner(Y = y_tune_STAARO, X = STAARO_Combined_Tune, family = gaussian(),
                        # For a real analysis we would use V = 10.
                        # V = 3,
                        SL.library = SL.library,cvControl = list(V = 20))

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

prs_columns <- c(which(str_detect(colnames(pheno_tuning_STAARO),"GeneCentric_Coding_")),which(str_detect(colnames(pheno_tuning_STAARO),"GeneCentric_Noncoding_")),which(str_detect(colnames(pheno_tuning_STAARO),"SlidingWindow_")),which(str_detect(colnames(pheno_tuning_STAARO),"prs")))

r2_tune <- vector()
for(j in 1:length(prs_columns)){
  r2_tune[j] <- summary(lm(as.formula(paste0("y_tune ~",colnames(pheno_tuning_STAARO)[prs_columns[j]])),data = pheno_tuning_STAARO))$r.squared
}
prs_best_tune <- data.frame(IID = pheno_tuning_STAARO$IID,prs = pheno_tuning_STAARO[,colnames(pheno_tuning_STAARO)[prs_columns[which.max(r2_tune)]]])

prs_best_validation <- data.frame(IID = pheno_vad_STAARO$IID,prs = pheno_vad_STAARO[,colnames(pheno_tuning_STAARO)[prs_columns[which.max(r2_tune)]]])

write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/Best_All_STAARO_Tune_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/Best_All_STAARO_Validation_All",i,".txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

y_vad_STAARO_EUR <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]]
y_vad_STAARO_NonEUR <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"]]
y_vad_STAARO_UNK <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"]]
y_vad_STAARO_SAS <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]]
y_vad_STAARO_MIX <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"]]
y_vad_STAARO_AFR <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]]
y_vad_STAARO_EAS <- y_vad_STAARO[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]]

prs_best_validation_EUR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
prs_best_validation_NonEur <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
prs_best_validation_UNK <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
prs_best_validation_SAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
prs_best_validation_MIX <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
prs_best_validation_AFR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
prs_best_validation_EAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]


## 
model <- lm(y_vad_STAARO_EUR~prs_best_validation_EUR$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_EUR, x = prs_best_validation_EUR$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_EUR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_EUR",i,".RData"))

## 
model <- lm(y_vad_STAARO_NonEUR~prs_best_validation_NonEur$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_NonEUR, x = prs_best_validation_NonEur$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_NonEUR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_NonEUR",i,".RData"))

## 
model <- lm(y_vad_STAARO_UNK~prs_best_validation_UNK$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_UNK, x = prs_best_validation_UNK$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_UNK",i,".RData"))

## 
model <- lm(y_vad_STAARO_UNK~prs_best_validation_UNK$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_UNK, x = prs_best_validation_UNK$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_UNK",i,".RData"))

## 
model <- lm(y_vad_STAARO_AFR~prs_best_validation_AFR$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_AFR, x = prs_best_validation_AFR$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_AFR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_AFR",i,".RData"))

## 
model <- lm(y_vad_STAARO_EAS~prs_best_validation_EAS$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_EAS, x = prs_best_validation_EAS$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_EAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_EAS",i,".RData"))

## 
model <- lm(y_vad_STAARO_SAS~prs_best_validation_SAS$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_SAS, x = prs_best_validation_SAS$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_SAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_SAS",i,".RData"))

## 
model <- lm(y_vad_STAARO_MIX~prs_best_validation_MIX$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_STAARO_MIX, x = prs_best_validation_MIX$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_STAARO_MIX",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_MIX",i,".RData"))






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

prs_columns <- c(which(str_detect(colnames(pheno_tuning_Burden),"GeneCentric_Coding_")),which(str_detect(colnames(pheno_tuning_Burden),"GeneCentric_Noncoding_")),which(str_detect(colnames(pheno_tuning_Burden),"SlidingWindow_")),which(str_detect(colnames(pheno_tuning_Burden),"prs")))

r2_tune <- vector()
for(j in 1:length(prs_columns)){
  r2_tune[j] <- summary(lm(as.formula(paste0("y_tune ~",colnames(pheno_tuning_Burden)[prs_columns[j]])),data = pheno_tuning_Burden))$r.squared
}
prs_best_tune <- data.frame(IID = pheno_tuning_Burden$IID,prs = pheno_tuning_Burden[,colnames(pheno_tuning_Burden)[prs_columns[which.max(r2_tune)]]])

prs_best_validation <- data.frame(IID = pheno_vad_Burden$IID,prs = pheno_vad_Burden[,colnames(pheno_tuning_Burden)[prs_columns[which.max(r2_tune)]]])

write.table(prs_best_tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/Best_All_Burden_Tune_All",i,".txt"),sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/Best_All_Burden_Validation_All",i,".txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

y_vad_Burden_EUR <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]]
y_vad_Burden_NonEUR <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"]]
y_vad_Burden_UNK <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"]]
y_vad_Burden_SAS <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]]
y_vad_Burden_MIX <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"]]
y_vad_Burden_AFR <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]]
y_vad_Burden_EAS <- y_vad_Burden[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]]

prs_best_validation_EUR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
prs_best_validation_NonEur <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
prs_best_validation_UNK <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
prs_best_validation_SAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
prs_best_validation_MIX <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
prs_best_validation_AFR <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
prs_best_validation_EAS <- prs_best_validation[prs_best_validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

## 
model <- lm(y_vad_Burden_EUR~prs_best_validation_EUR$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_EUR, x = prs_best_validation_EUR$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_EUR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_EUR",i,".RData"))

## 
model <- lm(y_vad_Burden_NonEUR~prs_best_validation_NonEur$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_NonEUR, x = prs_best_validation_NonEur$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_NonEUR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_NonEUR",i,".RData"))

## 
model <- lm(y_vad_Burden_UNK~prs_best_validation_UNK$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_UNK, x = prs_best_validation_UNK$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_UNK",i,".RData"))

## 
model <- lm(y_vad_Burden_UNK~prs_best_validation_UNK$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_UNK, x = prs_best_validation_UNK$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_UNK",i,".RData"))

## 
model <- lm(y_vad_Burden_AFR~prs_best_validation_AFR$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_AFR, x = prs_best_validation_AFR$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_AFR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_AFR",i,".RData"))

## 
model <- lm(y_vad_Burden_EAS~prs_best_validation_EAS$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_EAS, x = prs_best_validation_EAS$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_EAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_EAS",i,".RData"))

## 
model <- lm(y_vad_Burden_SAS~prs_best_validation_SAS$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_SAS, x = prs_best_validation_SAS$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_SAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_SAS",i,".RData"))

## 
model <- lm(y_vad_Burden_MIX~prs_best_validation_MIX$prs)
r2 <- summary(model)$r.square

data <- data.frame(y = y_vad_Burden_MIX, x = prs_best_validation_MIX$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Burden_MIX",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_MIX",i,".RData"))
  
