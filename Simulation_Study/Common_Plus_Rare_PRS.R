rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

## STAARO

pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"))

rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_STAARO_Tune_All",i,".txt"))

pheno_tune <- left_join(pheno_tune,common_prs,by = "IID")
colnames(pheno_tune) <- c(colnames(pheno_tune)[1:2],"CV_PRS")
pheno_tune <- left_join(pheno_tune,rarevariant_prs,by = "IID")
colnames(pheno_tune) <- c(colnames(pheno_tune)[1:3],"RV_PRS")

PRSs_Tune <- pheno_tune[,c(3:4)]

model.null <- lm(Y~1,data=pheno_tune)

PRSs_Tune$Residuals <- model.null$residuals

model.best <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))

rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_STAARO_Validation_All",i,".txt"))

pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
colnames(pheno_vad) <- c(colnames(pheno_vad)[1:2],"CV_PRS")
pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
colnames(pheno_vad) <- c(colnames(pheno_vad)[1:3],"RV_PRS")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]


##
PRSs_Validation <- pheno_vad_EUR[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_EUR)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_EUR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_EUR",i,".RData"))

##
PRSs_Validation <- pheno_vad_NonEur[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_NonEur)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_NonEur",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_NonEur",i,".RData"))

##
PRSs_Validation <- pheno_vad_UNK[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_UNK)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_UNK",i,".RData"))

##
PRSs_Validation <- pheno_vad_AFR[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_AFR)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_AFR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_AFR",i,".RData"))

##
PRSs_Validation <- pheno_vad_SAS[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_SAS)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_SAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_SAS",i,".RData"))

##
PRSs_Validation <- pheno_vad_MIX[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_MIX)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_MIX",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_MIX",i,".RData"))

##
PRSs_Validation <- pheno_vad_EAS[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_EAS)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_STAARO_EAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_EAS",i,".RData"))




## Burden

pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"))

rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_Burden_Tune_All",i,".txt"))

pheno_tune <- left_join(pheno_tune,common_prs,by = "IID")
colnames(pheno_tune) <- c(colnames(pheno_tune)[1:2],"CV_PRS")
pheno_tune <- left_join(pheno_tune,rarevariant_prs,by = "IID")
colnames(pheno_tune) <- c(colnames(pheno_tune)[1:3],"RV_PRS")

PRSs_Tune <- pheno_tune[,c(3,4)]

model.null <- lm(Y~1,data=pheno_tune)

PRSs_Tune$Residuals <- model.null$residuals

model.best <- lm(Residuals~CV_PRS + RV_PRS,data = PRSs_Tune)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))

rarevariant_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_All_Burden_Validation_All",i,".txt"))

pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
colnames(pheno_vad) <- c(colnames(pheno_vad)[1:2],"CV_PRS")
pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
colnames(pheno_vad) <- c(colnames(pheno_vad)[1:3],"RV_PRS")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

##
PRSs_Validation <- pheno_vad_EUR[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_EUR)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_EUR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_EUR",i,".RData"))

##
PRSs_Validation <- pheno_vad_NonEur[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_NonEur)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_NonEur",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_NonEur",i,".RData"))

##
PRSs_Validation <- pheno_vad_UNK[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_UNK)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_UNK",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_UNK",i,".RData"))

##
PRSs_Validation <- pheno_vad_AFR[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_AFR)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_AFR",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_AFR",i,".RData"))

##
PRSs_Validation <- pheno_vad_SAS[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_SAS)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_SAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_SAS",i,".RData"))

##
PRSs_Validation <- pheno_vad_MIX[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_MIX)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_MIX",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_MIX",i,".RData"))

##
PRSs_Validation <- pheno_vad_EAS[,c(3,4)]
model.null <- lm(Y~1,data=pheno_vad_EAS)
PRSs_Validation$Residuals <- model.null$residuals
predicted_prs <- predict(model.best,PRSs_Validation)
r2 <- summary(lm(PRSs_Validation$Residuals~predicted_prs))$r.square

data <- data.frame(y = PRSs_Validation$Residuals, x = predicted_prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "CV_plus_RV_Burden_EAS",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_EAS",i,".RData"))

  
