rm(list = ls())
library(dplyr)

STAARO_GeneCentric_Coding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/STAARO_GeneCentric_Coding_Tune_PRS.csv")
STAARO_GeneCentric_Coding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/STAARO_GeneCentric_Coding_Validation_PRS.csv")

STAARO_GeneCentric_Noncoding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/STAARO_GeneCentric_Noncoding_Tune_PRS.csv")
STAARO_GeneCentric_Noncoding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/STAARO_GeneCentric_Noncoding_Validation_PRS.csv")

STAARO_SlidingWindow_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/STAARO_SlidingWindow_Tune_PRS.csv")
STAARO_SlidingWindow_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/STAARO_SlidingWindow_Validation_PRS.csv")


Burden_GeneCentric_Coding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Burden_GeneCentric_Coding_Tune_PRS.csv")
Burden_GeneCentric_Coding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Burden_GeneCentric_Coding_Validation_PRS.csv")

Burden_GeneCentric_Noncoding_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Burden_GeneCentric_Noncoding_Tune_PRS.csv")
Burden_GeneCentric_Noncoding_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Burden_GeneCentric_Noncoding_Validation_PRS.csv")

Burden_SlidingWindow_Tune_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Burden_SlidingWindow_Tune_PRS.csv")
Burden_SlidingWindow_Validation_PRS <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Burden_SlidingWindow_Validation_PRS.csv")


## Pull in Phenotypes/Covariates 

pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
colnames(pheno_tuning) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

STAARO_GeneCentric_Coding_Tune_PRS <- left_join(pheno_tuning,STAARO_GeneCentric_Coding_Tune_PRS,by = "IID")
STAARO_GeneCentric_Noncoding_Tune_PRS <- left_join(pheno_tuning,STAARO_GeneCentric_Noncoding_Tune_PRS,by = "IID")
STAARO_SlidingWindow_Tune_PRS <- left_join(pheno_tuning,STAARO_SlidingWindow_Tune_PRS,by = "IID")

Burden_GeneCentric_Coding_Tune_PRS <- left_join(pheno_tuning,Burden_GeneCentric_Coding_Tune_PRS,by = "IID")
Burden_GeneCentric_Noncoding_Tune_PRS <- left_join(pheno_tuning,Burden_GeneCentric_Noncoding_Tune_PRS,by = "IID")
Burden_SlidingWindow_Tune_PRS <- left_join(pheno_tuning,Burden_SlidingWindow_Tune_PRS,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

STAARO_GeneCentric_Coding_Validation_PRS <- left_join(pheno_vad,STAARO_GeneCentric_Coding_Validation_PRS,by = "IID")
STAARO_GeneCentric_Noncoding_Validation_PRS <- left_join(pheno_vad,STAARO_GeneCentric_Noncoding_Validation_PRS,by = "IID")
STAARO_SlidingWindow_Validation_PRS <- left_join(pheno_vad,STAARO_SlidingWindow_Validation_PRS,by = "IID")

Burden_GeneCentric_Coding_Validation_PRS <- left_join(pheno_vad,Burden_GeneCentric_Coding_Validation_PRS,by = "IID")
Burden_GeneCentric_Noncoding_Validation_PRS <- left_join(pheno_vad,Burden_GeneCentric_Noncoding_Validation_PRS,by = "IID")
Burden_SlidingWindow_Validation_PRS <- left_join(pheno_vad,Burden_SlidingWindow_Validation_PRS,by = "IID")

############### R2's
arrayid <- as.numeric(commandArgs(TRUE)[1])


if(arrayid == 1){
  ##### Gene Centric Coding: STAARO
  r2_tun_GeneCentric_Coding_STAARO  <- rep(0,15)
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=STAARO_GeneCentric_Coding_Tune_PRS)
  for(k in 1:15){
    prs <- STAARO_GeneCentric_Coding_Tune_PRS[,paste0("PRS_Threshold_",k)]
    model.prs <- lm(model.null$residual~prs,data=STAARO_GeneCentric_Coding_Tune_PRS)
    r2_tun_GeneCentric_Coding_STAARO[k] <- summary(model.prs)$r.square
  }
  
  model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=STAARO_GeneCentric_Coding_Validation_PRS)
  prs <- STAARO_GeneCentric_Coding_Validation_PRS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_STAARO))]
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
  r2.result_GeneCentric_Coding_STAARO <- data.frame(method = "GeneCentric_Coding_STAARO",
                                                    r2 = r2_GeneCentric_Coding_STAARO,
                                                    r2_low = ci_result$percent[4],
                                                    r2_high = ci_result$percent[5]
  )
  
  save(r2.result_GeneCentric_Coding_STAARO, file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/GeneCentric_Coding_STAARO_result.RData") 
}else if(arrayid == 2){
  ##### Gene Centric Coding: Burden
  r2_tun_GeneCentric_Coding_Burden <- rep(0,15)
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=Burden_GeneCentric_Coding_Tune_PRS)
  for(k in 1:15){
    prs <- Burden_GeneCentric_Coding_Tune_PRS[,paste0("PRS_Threshold_",k)]
    model.prs <- lm(model.null$residual~prs,data=Burden_GeneCentric_Coding_Tune_PRS)
    r2_tun_GeneCentric_Coding_Burden[k] <- summary(model.prs)$r.square
  }
  
  model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=Burden_GeneCentric_Coding_Validation_PRS)
  prs <- Burden_GeneCentric_Coding_Validation_PRS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Coding_Burden))]
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
  r2.result_GeneCentric_Coding_Burden <- data.frame(method = "GeneCentric_Coding_Burden",
                                                    r2 = r2_GeneCentric_Coding_Burden,
                                                    r2_low = ci_result$percent[4],
                                                    r2_high = ci_result$percent[5]
  )
  
  save(r2.result_GeneCentric_Coding_Burden, file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/GeneCentric_Coding_Burden_result.RData")
}else if(arrayid == 3){
  ##### Gene Centric Noncoding: STAARO
  r2_tun_GeneCentric_Noncoding_STAARO  <- rep(0,15)
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=STAARO_GeneCentric_Noncoding_Tune_PRS)
  for(k in 1:15){
    prs <- STAARO_GeneCentric_Noncoding_Tune_PRS[,paste0("PRS_Threshold_",k)]
    model.prs <- lm(model.null$residual~prs,data=STAARO_GeneCentric_Noncoding_Tune_PRS)
    r2_tun_GeneCentric_Noncoding_STAARO[k] <- summary(model.prs)$r.square
  }
  
  model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=STAARO_GeneCentric_Noncoding_Validation_PRS)
  prs <- STAARO_GeneCentric_Noncoding_Validation_PRS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_STAARO))]
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
  r2.result_GeneCentric_Noncoding_STAARO <- data.frame(method = "GeneCentric_Noncoding_STAARO",
                                                       r2 = r2_GeneCentric_Noncoding_STAARO,
                                                       r2_low = ci_result$percent[4],
                                                       r2_high = ci_result$percent[5]
  )
  
  save(r2.result_GeneCentric_Noncoding_STAARO, file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/GeneCentric_Noncoding_STAARO_result.RData")
}else if(arrayid == 4){
  ##### Gene Centric Noncoding: Burden
  r2_tun_GeneCentric_Noncoding_Burden <- rep(0,15)
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=Burden_GeneCentric_Noncoding_Tune_PRS)
  for(k in 1:15){
    prs <- Burden_GeneCentric_Noncoding_Tune_PRS[,paste0("PRS_Threshold_",k)]
    model.prs <- lm(model.null$residual~prs,data=Burden_GeneCentric_Noncoding_Tune_PRS)
    r2_tun_GeneCentric_Noncoding_Burden[k] <- summary(model.prs)$r.square
  }
  
  model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=Burden_GeneCentric_Noncoding_Validation_PRS)
  prs <- Burden_GeneCentric_Noncoding_Validation_PRS[,paste0("PRS_Threshold_",which.max(r2_tun_GeneCentric_Noncoding_Burden))]
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
  r2.result_GeneCentric_Noncoding_Burden <- data.frame(method = "GeneCentric_Noncoding_Burden",
                                                       r2 = r2_GeneCentric_Noncoding_Burden,
                                                       r2_low = ci_result$percent[4],
                                                       r2_high = ci_result$percent[5]
  )
  
  save(r2.result_GeneCentric_Noncoding_Burden, file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/GeneCentric_Noncoding_Burden_result.RData")
}else if(arrayid == 5){
  ##### Gene Centric Sliding Window: STAARO
  r2_tun_SlidingWindow_STAARO  <- rep(0,15)
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=STAARO_SlidingWindow_Tune_PRS)
  for(k in 1:15){
    prs <- STAARO_SlidingWindow_Tune_PRS[,paste0("PRS_Threshold_",k)]
    model.prs <- lm(model.null$residual~prs,data=STAARO_SlidingWindow_Tune_PRS)
    r2_tun_SlidingWindow_STAARO[k] <- summary(model.prs)$r.square
  }
  
  model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=STAARO_SlidingWindow_Validation_PRS)
  prs <- STAARO_SlidingWindow_Validation_PRS[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_STAARO))]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2_SlidingWindow_STAARO <- summary(model.vad.prs)$r.square
  
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
  r2.result_SlidingWindow_STAARO <- data.frame(method = "SlidingWindow_STAARO",
                                               r2 = r2_SlidingWindow_STAARO,
                                               r2_low = ci_result$percent[4],
                                               r2_high = ci_result$percent[5]
  )
  
  save(r2.result_SlidingWindow_STAARO, file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/SlidingWindow_STAARO_result.RData")
}else{
  ##### Gene Centric Sliding Window: Burden
  r2_tun_SlidingWindow_Burden <- rep(0,15)
  model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=Burden_SlidingWindow_Tune_PRS)
  for(k in 1:15){
    prs <- Burden_SlidingWindow_Tune_PRS[,paste0("PRS_Threshold_",k)]
    model.prs <- lm(model.null$residual~prs,data=Burden_SlidingWindow_Tune_PRS)
    r2_tun_SlidingWindow_Burden[k] <- summary(model.prs)$r.square
  }
  
  model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=Burden_SlidingWindow_Validation_PRS)
  prs <- Burden_SlidingWindow_Validation_PRS[,paste0("PRS_Threshold_",which.max(r2_tun_SlidingWindow_Burden))]
  model.vad.prs <- lm(model.vad.null$residual~prs)
  r2_SlidingWindow_Burden <- summary(model.vad.prs)$r.square
  
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
  r2.result_SlidingWindow_Burden <- data.frame(method = "SlidingWindow_Burden",
                                               r2 = r2_SlidingWindow_Burden,
                                               r2_low = ci_result$percent[4],
                                               r2_high = ci_result$percent[5]
  )
  
  save(r2.result_SlidingWindow_Burden, file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/SlidingWindow_Burden_result.RData")
}






