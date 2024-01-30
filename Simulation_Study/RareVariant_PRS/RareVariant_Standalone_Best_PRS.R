rm(list = ls())
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

# load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/tune_prs_mat",i,".RData"))
# prs_mat <- data.frame(IID = ids_tune,PRS = prs_mat)
# STAARO_SlidingWindow_Tune_PRS <- prs_mat[,c(1,2:16)]
# colnames(STAARO_SlidingWindow_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
# Burden_SlidingWindow_Tune_PRS <- prs_mat[,c(1,17:31)]
# colnames(Burden_SlidingWindow_Tune_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
# 
# load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/validation_prs_mat",i,".RData"))
# prs_mat <- data.frame(IID = ids_validation,PRS = prs_mat)
# STAARO_SlidingWindow_Validation_PRS <- prs_mat[,c(1,2:16)]
# colnames(STAARO_SlidingWindow_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))
# Burden_SlidingWindow_Validation_PRS <- prs_mat[,c(1,17:31)]
# colnames(Burden_SlidingWindow_Validation_PRS) <- c("IID",paste0("PRS_Threshold_",1:15))

## Pull in Phenotypes/Covariates 

pheno_tuning <- Y_tune[[i]]
colnames(pheno_tuning) <- c("IID","Y")
STAARO_GeneCentric_Coding_Tune_PRS <- left_join(pheno_tuning,STAARO_GeneCentric_Coding_Tune_PRS,by = "IID")
STAARO_GeneCentric_Noncoding_Tune_PRS <- left_join(pheno_tuning,STAARO_GeneCentric_Noncoding_Tune_PRS,by = "IID")
# STAARO_SlidingWindow_Tune_PRS <- left_join(pheno_tuning,STAARO_SlidingWindow_Tune_PRS,by = "IID")

Burden_GeneCentric_Coding_Tune_PRS <- left_join(pheno_tuning,Burden_GeneCentric_Coding_Tune_PRS,by = "IID")
Burden_GeneCentric_Noncoding_Tune_PRS <- left_join(pheno_tuning,Burden_GeneCentric_Noncoding_Tune_PRS,by = "IID")
# Burden_SlidingWindow_Tune_PRS <- left_join(pheno_tuning,Burden_SlidingWindow_Tune_PRS,by = "IID")

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

STAARO_GeneCentric_Coding_Validation_PRS <- left_join(pheno_vad,STAARO_GeneCentric_Coding_Validation_PRS,by = "IID")
STAARO_GeneCentric_Coding_Validation_PRS_EUR <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
STAARO_GeneCentric_Coding_Validation_PRS_NonEur <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
STAARO_GeneCentric_Coding_Validation_PRS_UNK <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
STAARO_GeneCentric_Coding_Validation_PRS_SAS <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
STAARO_GeneCentric_Coding_Validation_PRS_MIX <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
STAARO_GeneCentric_Coding_Validation_PRS_AFR <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
STAARO_GeneCentric_Coding_Validation_PRS_EAS <- STAARO_GeneCentric_Coding_Validation_PRS[STAARO_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

STAARO_GeneCentric_Noncoding_Validation_PRS <- left_join(pheno_vad,STAARO_GeneCentric_Noncoding_Validation_PRS,by = "IID")
STAARO_GeneCentric_Noncoding_Validation_PRS_EUR <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
STAARO_GeneCentric_Noncoding_Validation_PRS_NonEur <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
STAARO_GeneCentric_Noncoding_Validation_PRS_UNK <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
STAARO_GeneCentric_Noncoding_Validation_PRS_SAS <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
STAARO_GeneCentric_Noncoding_Validation_PRS_MIX <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
STAARO_GeneCentric_Noncoding_Validation_PRS_AFR <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
STAARO_GeneCentric_Noncoding_Validation_PRS_EAS <- STAARO_GeneCentric_Noncoding_Validation_PRS[STAARO_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

# STAARO_SlidingWindow_Validation_PRS <- left_join(pheno_vad,STAARO_SlidingWindow_Validation_PRS,by = "IID")
# STAARO_SlidingWindow_Validation_PRS_EUR <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
# STAARO_SlidingWindow_Validation_PRS_NonEur <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
# STAARO_SlidingWindow_Validation_PRS_UNK <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
# STAARO_SlidingWindow_Validation_PRS_SAS <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
# STAARO_SlidingWindow_Validation_PRS_MIX <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
# STAARO_SlidingWindow_Validation_PRS_AFR <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
# STAARO_SlidingWindow_Validation_PRS_EAS <- STAARO_SlidingWindow_Validation_PRS[STAARO_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

Burden_GeneCentric_Coding_Validation_PRS <- left_join(pheno_vad,Burden_GeneCentric_Coding_Validation_PRS,by = "IID")
Burden_GeneCentric_Coding_Validation_PRS_EUR <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
Burden_GeneCentric_Coding_Validation_PRS_NonEur <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
Burden_GeneCentric_Coding_Validation_PRS_UNK <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
Burden_GeneCentric_Coding_Validation_PRS_SAS <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
Burden_GeneCentric_Coding_Validation_PRS_MIX <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
Burden_GeneCentric_Coding_Validation_PRS_AFR <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
Burden_GeneCentric_Coding_Validation_PRS_EAS <- Burden_GeneCentric_Coding_Validation_PRS[Burden_GeneCentric_Coding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

Burden_GeneCentric_Noncoding_Validation_PRS <- left_join(pheno_vad,Burden_GeneCentric_Noncoding_Validation_PRS,by = "IID")
Burden_GeneCentric_Noncoding_Validation_PRS_EUR <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
Burden_GeneCentric_Noncoding_Validation_PRS_NonEur <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
Burden_GeneCentric_Noncoding_Validation_PRS_UNK <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
Burden_GeneCentric_Noncoding_Validation_PRS_SAS <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
Burden_GeneCentric_Noncoding_Validation_PRS_MIX <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
Burden_GeneCentric_Noncoding_Validation_PRS_AFR <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
Burden_GeneCentric_Noncoding_Validation_PRS_EAS <- Burden_GeneCentric_Noncoding_Validation_PRS[Burden_GeneCentric_Noncoding_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

# Burden_SlidingWindow_Validation_PRS <- left_join(pheno_vad,Burden_SlidingWindow_Validation_PRS,by = "IID")
# Burden_SlidingWindow_Validation_PRS_EUR <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
# Burden_SlidingWindow_Validation_PRS_NonEur <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
# Burden_SlidingWindow_Validation_PRS_UNK <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
# Burden_SlidingWindow_Validation_PRS_SAS <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
# Burden_SlidingWindow_Validation_PRS_MIX <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
# Burden_SlidingWindow_Validation_PRS_AFR <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
# Burden_SlidingWindow_Validation_PRS_EAS <- Burden_SlidingWindow_Validation_PRS[Burden_SlidingWindow_Validation_PRS$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

############### R2's

##### Gene Centric Coding: STAARO
r2_tun_GeneCentric_Coding_STAARO  <- rep(0,15)
model.null <- lm(Y~1,data=STAARO_GeneCentric_Coding_Tune_PRS)
for(k in 1:15){
  prs <- STAARO_GeneCentric_Coding_Tune_PRS[,paste0("PRS_Threshold_",k)]
  model.prs <- lm(model.null$residual~prs,data=STAARO_GeneCentric_Coding_Tune_PRS)
  r2_tun_GeneCentric_Coding_STAARO[k] <- summary(model.prs)$r.square
}


####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_EUR)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_EUR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EUR",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_EUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_NonEur)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_NonEur[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_NonEUR",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_NonEUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_UNK)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_UNK[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_UNK",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_UNK",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_AFR)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_AFR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_AFR",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_AFR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_SAS)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_SAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_SAS",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_SAS",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_MIX)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_MIX[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_MIX",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_MIX",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Coding_Validation_PRS_EAS)
prs <- STAARO_GeneCentric_Coding_Validation_PRS_EAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO_EAS",
                                                  r2 = r2_GeneCentric_Coding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_EAS",i,".RData"))

##### Gene Centric Coding: Burden
r2_tun_GeneCentric_Coding_Burden <- rep(0,15)
model.null <- lm(Y~1,data=Burden_GeneCentric_Coding_Tune_PRS)
for(k in 1:15){
  prs <- Burden_GeneCentric_Coding_Tune_PRS[,paste0("PRS_Threshold_",k)]
  model.prs <- lm(model.null$residual~prs,data=Burden_GeneCentric_Coding_Tune_PRS)
  r2_tun_GeneCentric_Coding_Burden[k] <- summary(model.prs)$r.square
}

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_EUR)
prs <- Burden_GeneCentric_Coding_Validation_PRS_EUR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EUR",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_EUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_NonEur)
prs <- Burden_GeneCentric_Coding_Validation_PRS_NonEur[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_NonEUR",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_NonEUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_UNK)
prs <- Burden_GeneCentric_Coding_Validation_PRS_UNK[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_UNK",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_UNK",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_AFR)
prs <- Burden_GeneCentric_Coding_Validation_PRS_AFR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_AFR",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_AFR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_SAS)
prs <- Burden_GeneCentric_Coding_Validation_PRS_SAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_SAS",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_SAS",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_MIX)
prs <- Burden_GeneCentric_Coding_Validation_PRS_MIX[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_MIX",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_MIX",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Coding_Validation_PRS_EAS)
prs <- Burden_GeneCentric_Coding_Validation_PRS_EAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Coding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden_EAS",
                                                  r2 = r2_GeneCentric_Coding_Burden,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Coding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_EAS",i,".RData"))


##### Gene Centric Noncoding: STAARO
r2_tun_GeneCentric_Noncoding_STAARO  <- rep(0,15)
model.null <- lm(Y~1,data=STAARO_GeneCentric_Noncoding_Tune_PRS)
for(k in 1:15){
  prs <- STAARO_GeneCentric_Noncoding_Tune_PRS[,paste0("PRS_Threshold_",k)]
  model.prs <- lm(model.null$residual~prs,data=STAARO_GeneCentric_Noncoding_Tune_PRS)
  r2_tun_GeneCentric_Noncoding_STAARO[k] <- summary(model.prs)$r.square
}

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_EUR)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_EUR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EUR",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_EUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_NonEur)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_NonEur[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_NonEUR",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_NonEUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_UNK)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_UNK[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_UNK",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_UNK",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_AFR)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_AFR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_AFR",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_AFR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_SAS)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_SAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_SAS",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_SAS",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_MIX)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_MIX[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_MIX",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_MIX",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=STAARO_GeneCentric_Noncoding_Validation_PRS_EAS)
prs <- STAARO_GeneCentric_Noncoding_Validation_PRS_EAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_STAARO <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO_EAS",
                                                  r2 = r2_GeneCentric_Noncoding_STAARO,
                                                  r2_low = ci_result$percent[4],
                                                  r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_EAS",i,".RData"))

##### Gene Centric Noncoding: Burden
r2_tun_GeneCentric_Noncoding_Burden <- rep(0,15)
model.null <- lm(Y~1,data=Burden_GeneCentric_Noncoding_Tune_PRS)
for(k in 1:15){
  prs <- Burden_GeneCentric_Noncoding_Tune_PRS[,paste0("PRS_Threshold_",k)]
  model.prs <- lm(model.null$residual~prs,data=Burden_GeneCentric_Noncoding_Tune_PRS)
  r2_tun_GeneCentric_Noncoding_Burden[k] <- summary(model.prs)$r.square
}

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_EUR)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_EUR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EUR",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_EUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_NonEur)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_NonEur[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_NonEUR",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_NonEUR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_UNK)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_UNK[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_UNK",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_UNK",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_AFR)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_AFR[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_AFR",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_AFR",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_SAS)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_SAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_SAS",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_SAS",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_MIX)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_MIX[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_MIX",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_MIX",i,".RData"))

####
model.vad.null  <-  lm(Y~1,data=Burden_GeneCentric_Noncoding_Validation_PRS_EAS)
prs <- Burden_GeneCentric_Noncoding_Validation_PRS_EAS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2_GeneCentric_Noncoding_Burden <- summary(model.vad.prs)$r.square

data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden_EAS",
                                                     r2 = r2_GeneCentric_Noncoding_Burden,
                                                     r2_low = ci_result$percent[4],
                                                     r2_high = ci_result$percent[5]
)

save(r2.result_GeneCentric_Noncoding_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_EAS",i,".RData"))

# ##### Gene Centric Sliding Window: STAARO
# r2_tun_SlidingWindow_STAARO  <- rep(0,15)
# model.null <- lm(Y~1,data=STAARO_SlidingWindow_Tune_PRS)
# for(k in 1:15){
#   prs <- STAARO_SlidingWindow_Tune_PRS[,paste0("PRS_Threshold_",k)]
#   model.prs <- lm(model.null$residual~prs,data=STAARO_SlidingWindow_Tune_PRS)
#   r2_tun_SlidingWindow_STAARO[k] <- summary(model.prs)$r.square
# }
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_EUR)
# prs <- STAARO_SlidingWindow_Validation_PRS_EUR[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_EUR",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_EUR",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_NonEur)
# prs <- STAARO_SlidingWindow_Validation_PRS_NonEur[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_NonEUR",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_NonEUR",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_UNK)
# prs <- STAARO_SlidingWindow_Validation_PRS_UNK[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_UNK",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_UNK",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_AFR)
# prs <- STAARO_SlidingWindow_Validation_PRS_AFR[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_AFR",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_AFR",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_SAS)
# prs <- STAARO_SlidingWindow_Validation_PRS_SAS[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_SAS",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_SAS",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_MIX)
# prs <- STAARO_SlidingWindow_Validation_PRS_MIX[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_MIX",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_MIX",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=STAARO_SlidingWindow_Validation_PRS_EAS)
# prs <- STAARO_SlidingWindow_Validation_PRS_EAS[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO_EAS",
#                                                      r2 = r2_SlidingWindow_STAARO,
#                                                      r2_low = ci_result$percent[4],
#                                                      r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_STAARO, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_EAS",i,".RData"))
# 
# ##### Gene Centric Sliding Window: Burden
# r2_tun_SlidingWindow_Burden <- rep(0,15)
# model.null <- lm(Y~1,data=Burden_SlidingWindow_Tune_PRS)
# for(k in 1:15){
#   prs <- Burden_SlidingWindow_Tune_PRS[,paste0("PRS_Threshold_",k)]
#   model.prs <- lm(model.null$residual~prs,data=Burden_SlidingWindow_Tune_PRS)
#   r2_tun_SlidingWindow_Burden[k] <- summary(model.prs)$r.square
# }
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_EUR)
# prs <- Burden_SlidingWindow_Validation_PRS_EUR[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_EUR",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_EUR",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_NonEur)
# prs <- Burden_SlidingWindow_Validation_PRS_NonEur[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_NonEUR",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_NonEUR",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_UNK)
# prs <- Burden_SlidingWindow_Validation_PRS_UNK[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_UNK",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_UNK",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_AFR)
# prs <- Burden_SlidingWindow_Validation_PRS_AFR[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_AFR",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_AFR",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_SAS)
# prs <- Burden_SlidingWindow_Validation_PRS_SAS[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_SAS",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_SAS",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_MIX)
# prs <- Burden_SlidingWindow_Validation_PRS_MIX[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_MIX",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_MIX",i,".RData"))
# 
# ####
# model.vad.null  <-  lm(Y~1,data=Burden_SlidingWindow_Validation_PRS_EAS)
# prs <- Burden_SlidingWindow_Validation_PRS_EAS[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
# model.vad.prs <- lm(model.vad.null$residual~prs)
# r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
# 
# data <- data.frame(y = model.vad.null$residual, x = prs)
# R2Boot <- function(data,indices){
#   boot_data <- data[indices, ]
#   model <- lm(y ~ x, data = boot_data)
#   result <- summary(model)$r.square
#   return(c(result))
# }
# library(boot)
# boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)
# 
# ci_result <- boot.ci(boot_r2, type = "perc")
# r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden_EAS",
#                                              r2 = r2_SlidingWindow_Burden,
#                                              r2_low = ci_result$percent[4],
#                                              r2_high = ci_result$percent[5]
# )
# 
# save(r2.result_SlidingWindow_Burden, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_EAS",i,".RData"))
# 
