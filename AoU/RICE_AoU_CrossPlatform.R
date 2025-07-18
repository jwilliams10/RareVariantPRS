rm(list = ls())

library(data.table)
library(dplyr)
library(boot)

trait <- c("Height","BMI","TC","HDL","LDL","logTG")[as.numeric(commandArgs(TRUE)[1])]

RICE_CV_PRS <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_RICECV_PRS.sscore"), header=FALSE, comment.char="#")
RICE_CV_PRS <- RICE_CV_PRS[,c(2,5)]
colnames(RICE_CV_PRS) <- c("IID","CV_PRS")

RICE_RV_PRS <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_RICERV_PRS.csv"))
colnames(RICE_RV_PRS) <- c("IID","RV_PRS")

RICE_PRS <- inner_join(RICE_CV_PRS,RICE_RV_PRS)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
ukb_pheno <- as.data.frame(ukb_pheno)
sampleids_all <- read.table("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega.fam", quote="\"", comment.char="")
unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")
ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids_all[,2],]
ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% unrels_nRandomSNPs_0[,2],]

ukb_pheno <- ukb_pheno[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
ukb_pheno$age2 <- ukb_pheno$age^2
ukb_pheno[c(14,16:26)] <- lapply(ukb_pheno[c(14,16:26)], function(x){c(scale(x))})

model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=ukb_pheno)
ukb_pheno$y_residualized <- NA
ukb_pheno$y_residualized[!is.na(ukb_pheno[,trait])] <- model.null$residual

ukb_pheno <- inner_join(ukb_pheno,RICE_PRS)

ukb_pheno_raw <- ukb_pheno[,c("IID","y_residualized","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","CV_PRS","RV_PRS")]
ukb_pheno_adjusted <- ukb_pheno_raw

for(i in c("CV_PRS","RV_PRS")){
  tmp <- data.frame(y = ukb_pheno_adjusted[,i],ukb_pheno_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
  mod <- lm(y~.,data = tmp)
  R <- mod$residuals
  tmp <- data.frame(y = R^2,ukb_pheno_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
  mod <- lm(y~.,data = tmp)
  y_hat <- predict(mod,tmp)
  if(sum(y_hat < 0) > 0){
    mod <- lm(y~1,data = tmp)
    y_hat <- predict(mod,tmp)
  }
  if(sum(sqrt(y_hat)) == 0){
    ukb_pheno_adjusted[,i] <- 0
  }else{
    ukb_pheno_adjusted[,i] <- R/sqrt(y_hat)
  }
}

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

ukb_pheno_raw_EUR <- ukb_pheno_raw[ukb_pheno_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
ukb_pheno_raw_SAS <- ukb_pheno_raw[ukb_pheno_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
ukb_pheno_raw_AMR <- ukb_pheno_raw[ukb_pheno_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
ukb_pheno_raw_AFR <- ukb_pheno_raw[ukb_pheno_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
ukb_pheno_raw_EAS <- ukb_pheno_raw[ukb_pheno_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

ukb_pheno_adjusted_EUR <- ukb_pheno_adjusted[ukb_pheno_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
ukb_pheno_adjusted_SAS <- ukb_pheno_adjusted[ukb_pheno_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
ukb_pheno_adjusted_AMR <- ukb_pheno_adjusted[ukb_pheno_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
ukb_pheno_adjusted_AFR <- ukb_pheno_adjusted[ukb_pheno_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
ukb_pheno_adjusted_EAS <- ukb_pheno_adjusted[ukb_pheno_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

ukb_pheno_raw_EUR$y_residualized <- scale(ukb_pheno_raw_EUR$y_residualized)
ukb_pheno_raw_SAS$y_residualized <- scale(ukb_pheno_raw_SAS$y_residualized)
ukb_pheno_raw_AMR$y_residualized <- scale(ukb_pheno_raw_AMR$y_residualized)
ukb_pheno_raw_AFR$y_residualized <- scale(ukb_pheno_raw_AFR$y_residualized)
ukb_pheno_raw_EAS$y_residualized <- scale(ukb_pheno_raw_EAS$y_residualized)

ukb_pheno_adjusted_EUR$y_residualized <- scale(ukb_pheno_adjusted_EUR$y_residualized)
ukb_pheno_adjusted_SAS$y_residualized <- scale(ukb_pheno_adjusted_SAS$y_residualized)
ukb_pheno_adjusted_AMR$y_residualized <- scale(ukb_pheno_adjusted_AMR$y_residualized)
ukb_pheno_adjusted_AFR$y_residualized <- scale(ukb_pheno_adjusted_AFR$y_residualized)
ukb_pheno_adjusted_EAS$y_residualized <- scale(ukb_pheno_adjusted_EAS$y_residualized)

ukb_pheno_raw_EUR[,"CV_PRS"] <- scale(ukb_pheno_raw_EUR[,"CV_PRS"])
ukb_pheno_raw_SAS[,"CV_PRS"] <- scale(ukb_pheno_raw_SAS[,"CV_PRS"])
ukb_pheno_raw_AMR[,"CV_PRS"] <- scale(ukb_pheno_raw_AMR[,"CV_PRS"])
ukb_pheno_raw_AFR[,"CV_PRS"] <- scale(ukb_pheno_raw_AFR[,"CV_PRS"])
ukb_pheno_raw_EAS[,"CV_PRS"] <- scale(ukb_pheno_raw_EAS[,"CV_PRS"])

ukb_pheno_raw_EUR[,"RV_PRS"] <- scale(ukb_pheno_raw_EUR[,"RV_PRS"])
ukb_pheno_raw_SAS[,"RV_PRS"] <- scale(ukb_pheno_raw_SAS[,"RV_PRS"])
ukb_pheno_raw_AMR[,"RV_PRS"] <- scale(ukb_pheno_raw_AMR[,"RV_PRS"])
ukb_pheno_raw_AFR[,"RV_PRS"] <- scale(ukb_pheno_raw_AFR[,"RV_PRS"])
ukb_pheno_raw_EAS[,"RV_PRS"] <- scale(ukb_pheno_raw_EAS[,"RV_PRS"])

Beta_CV_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = boot_data))[2]
  return(c(result))
}

Beta_RV_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = boot_data))[3]
  return(c(result))
}

R2_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = boot_data))$r.squared
  return(c(result))
}

beta_CV_validation_raw_EUR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_EUR))[2]
boot_beta <- boot(data = ukb_pheno_raw_EUR, statistic = Beta_CV_Boot, R = 10000)
beta_CV_raw_EUR_boot <- boot_beta$t
beta_CV_se_validation_raw_EUR <- sd(boot_beta$t)

beta_RV_validation_raw_EUR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_EUR))[3]
boot_beta <- boot(data = ukb_pheno_raw_EUR, statistic = Beta_RV_Boot, R = 10000)
beta_RV_raw_EUR_boot <- boot_beta$t
beta_RV_se_validation_raw_EUR <- sd(boot_beta$t)

R2_validation_raw_EUR <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_EUR))$r.squared
boot_R2 <- boot(data = ukb_pheno_raw_EUR, statistic = R2_Boot, R = 10000)
R2_raw_EUR_boot <- boot_R2$t
R2_se_validation_raw_EUR <- sd(boot_R2$t)

beta_CV_validation_raw_SAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_SAS))[2]
boot_beta <- boot(data = ukb_pheno_raw_SAS, statistic = Beta_CV_Boot, R = 10000)
beta_CV_raw_SAS_boot <- boot_beta$t
beta_CV_se_validation_raw_SAS <- sd(boot_beta$t)

beta_RV_validation_raw_SAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_SAS))[3]
boot_beta <- boot(data = ukb_pheno_raw_SAS, statistic = Beta_RV_Boot, R = 10000)
beta_RV_raw_SAS_boot <- boot_beta$t
beta_RV_se_validation_raw_SAS <- sd(boot_beta$t)

R2_validation_raw_SAS <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_SAS))$r.squared
boot_R2 <- boot(data = ukb_pheno_raw_SAS, statistic = R2_Boot, R = 10000)
R2_raw_SAS_boot <- boot_R2$t
R2_se_validation_raw_SAS <- sd(boot_R2$t)

beta_CV_validation_raw_AMR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_AMR))[2]
boot_beta <- boot(data = ukb_pheno_raw_AMR, statistic = Beta_CV_Boot, R = 10000)
beta_CV_raw_AMR_boot <- boot_beta$t
beta_CV_se_validation_raw_AMR <- sd(boot_beta$t)

beta_RV_validation_raw_AMR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_AMR))[3]
boot_beta <- boot(data = ukb_pheno_raw_AMR, statistic = Beta_RV_Boot, R = 10000)
beta_RV_raw_AMR_boot <- boot_beta$t
beta_RV_se_validation_raw_AMR <- sd(boot_beta$t)

R2_validation_raw_AMR <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_AMR))$r.squared
boot_R2 <- boot(data = ukb_pheno_raw_AMR, statistic = R2_Boot, R = 10000)
R2_raw_AMR_boot <- boot_R2$t
R2_se_validation_raw_AMR <- sd(boot_R2$t)

beta_CV_validation_raw_AFR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_AFR))[2]
boot_beta <- boot(data = ukb_pheno_raw_AFR, statistic = Beta_CV_Boot, R = 10000)
beta_CV_raw_AFR_boot <- boot_beta$t
beta_CV_se_validation_raw_AFR <- sd(boot_beta$t)

beta_RV_validation_raw_AFR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_AFR))[3]
boot_beta <- boot(data = ukb_pheno_raw_AFR, statistic = Beta_RV_Boot, R = 10000)
beta_RV_raw_AFR_boot <- boot_beta$t
beta_RV_se_validation_raw_AFR <- sd(boot_beta$t)

R2_validation_raw_AFR <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_AFR))$r.squared
boot_R2 <- boot(data = ukb_pheno_raw_AFR, statistic = R2_Boot, R = 10000)
R2_raw_AFR_boot <- boot_R2$t
R2_se_validation_raw_AFR <- sd(boot_R2$t)

beta_CV_validation_raw_EAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_EAS))[2]
boot_beta <- boot(data = ukb_pheno_raw_EAS, statistic = Beta_CV_Boot, R = 10000)
beta_CV_raw_EAS_boot <- boot_beta$t
beta_CV_se_validation_raw_EAS <- sd(boot_beta$t)

beta_RV_validation_raw_EAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_EAS))[3]
boot_beta <- boot(data = ukb_pheno_raw_EAS, statistic = Beta_RV_Boot, R = 10000)
beta_RV_raw_EAS_boot <- boot_beta$t
beta_RV_se_validation_raw_EAS <- sd(boot_beta$t)

R2_validation_raw_EAS <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_raw_EAS))$r.squared
boot_R2 <- boot(data = ukb_pheno_raw_EAS, statistic = R2_Boot, R = 10000)
R2_raw_EAS_boot <- boot_R2$t
R2_se_validation_raw_EAS <- sd(boot_R2$t)

beta_CV_validation_adjusted_EUR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_EUR))[2]
boot_beta <- boot(data = ukb_pheno_adjusted_EUR, statistic = Beta_CV_Boot, R = 10000)
beta_CV_adjusted_EUR_boot <- boot_beta$t
beta_CV_se_validation_adjusted_EUR <- sd(boot_beta$t)

beta_RV_validation_adjusted_EUR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_EUR))[3]
boot_beta <- boot(data = ukb_pheno_adjusted_EUR, statistic = Beta_RV_Boot, R = 10000)
beta_RV_adjusted_EUR_boot <- boot_beta$t
beta_RV_se_validation_adjusted_EUR <- sd(boot_beta$t)

R2_validation_adjusted_EUR <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_EUR))$r.squared
boot_R2 <- boot(data = ukb_pheno_adjusted_EUR, statistic = R2_Boot, R = 10000)
R2_adjusted_EUR_boot <- boot_R2$t
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)

beta_CV_validation_adjusted_SAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_SAS))[2]
boot_beta <- boot(data = ukb_pheno_adjusted_SAS, statistic = Beta_CV_Boot, R = 10000)
beta_CV_adjusted_SAS_boot <- boot_beta$t
beta_CV_se_validation_adjusted_SAS <- sd(boot_beta$t)

beta_RV_validation_adjusted_SAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_SAS))[3]
boot_beta <- boot(data = ukb_pheno_adjusted_SAS, statistic = Beta_RV_Boot, R = 10000)
beta_RV_adjusted_SAS_boot <- boot_beta$t
beta_RV_se_validation_adjusted_SAS <- sd(boot_beta$t)

R2_validation_adjusted_SAS <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_SAS))$r.squared
boot_R2 <- boot(data = ukb_pheno_adjusted_SAS, statistic = R2_Boot, R = 10000)
R2_adjusted_SAS_boot <- boot_R2$t
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)

beta_CV_validation_adjusted_AMR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_AMR))[2]
boot_beta <- boot(data = ukb_pheno_adjusted_AMR, statistic = Beta_CV_Boot, R = 10000)
beta_CV_adjusted_AMR_boot <- boot_beta$t
beta_CV_se_validation_adjusted_AMR <- sd(boot_beta$t)

beta_RV_validation_adjusted_AMR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_AMR))[3]
boot_beta <- boot(data = ukb_pheno_adjusted_AMR, statistic = Beta_RV_Boot, R = 10000)
beta_RV_adjusted_AMR_boot <- boot_beta$t
beta_RV_se_validation_adjusted_AMR <- sd(boot_beta$t)

R2_validation_adjusted_AMR <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_AMR))$r.squared
boot_R2 <- boot(data = ukb_pheno_adjusted_AMR, statistic = R2_Boot, R = 10000)
R2_adjusted_AMR_boot <- boot_R2$t
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)

beta_CV_validation_adjusted_AFR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_AFR))[2]
boot_beta <- boot(data = ukb_pheno_adjusted_AFR, statistic = Beta_CV_Boot, R = 10000)
beta_CV_adjusted_AFR_boot <- boot_beta$t
beta_CV_se_validation_adjusted_AFR <- sd(boot_beta$t)

beta_RV_validation_adjusted_AFR <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_AFR))[3]
boot_beta <- boot(data = ukb_pheno_adjusted_AFR, statistic = Beta_RV_Boot, R = 10000)
beta_RV_adjusted_AFR_boot <- boot_beta$t
beta_RV_se_validation_adjusted_AFR <- sd(boot_beta$t)

R2_validation_adjusted_AFR <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_AFR))$r.squared
boot_R2 <- boot(data = ukb_pheno_adjusted_AFR, statistic = R2_Boot, R = 10000)
R2_adjusted_AFR_boot <- boot_R2$t
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)

beta_CV_validation_adjusted_EAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_EAS))[2]
boot_beta <- boot(data = ukb_pheno_adjusted_EAS, statistic = Beta_CV_Boot, R = 10000)
beta_CV_adjusted_EAS_boot <- boot_beta$t
beta_CV_se_validation_adjusted_EAS <- sd(boot_beta$t)

beta_RV_validation_adjusted_EAS <- coef(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_EAS))[3]
boot_beta <- boot(data = ukb_pheno_adjusted_EAS, statistic = Beta_RV_Boot, R = 10000)
beta_RV_adjusted_EAS_boot <- boot_beta$t
beta_RV_se_validation_adjusted_EAS <- sd(boot_beta$t)

R2_validation_adjusted_EAS <- summary(lm(y_residualized~CV_PRS + RV_PRS,data = ukb_pheno_adjusted_EAS))$r.squared
boot_R2 <- boot(data = ukb_pheno_adjusted_EAS, statistic = R2_Boot, R = 10000)
R2_adjusted_EAS_boot <- boot_R2$t
R2_se_validation_adjusted_EAS <- sd(boot_R2$t)

CV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(beta_CV_validation_raw_EUR,beta_CV_validation_raw_SAS,beta_CV_validation_raw_AMR,beta_CV_validation_raw_AFR,beta_CV_validation_raw_EAS), 
                             beta_se_raw = c(beta_CV_se_validation_raw_EUR,beta_CV_se_validation_raw_SAS,beta_CV_se_validation_raw_AMR,beta_CV_se_validation_raw_AFR,beta_CV_se_validation_raw_EAS), 
                             R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                             R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                             beta_adjusted = c(beta_CV_validation_adjusted_EUR,beta_CV_validation_adjusted_SAS,beta_CV_validation_adjusted_AMR,beta_CV_validation_adjusted_AFR,beta_CV_validation_adjusted_EAS), 
                             beta_se_adjusted = c(beta_CV_se_validation_adjusted_EUR,beta_CV_se_validation_adjusted_SAS,beta_CV_se_validation_adjusted_AMR,beta_CV_se_validation_adjusted_AFR,beta_CV_se_validation_adjusted_EAS), 
                             R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                             R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))

CV_Boot_Results <- data.frame(trait = trait,beta_CV_raw_EUR_boot,R2_raw_EUR_boot,beta_CV_raw_SAS_boot,R2_raw_SAS_boot,
                              beta_CV_raw_AMR_boot,R2_raw_AMR_boot,beta_CV_raw_AFR_boot,R2_raw_AFR_boot,
                              beta_CV_raw_EAS_boot,R2_raw_EAS_boot,beta_CV_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                              beta_CV_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_CV_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                              beta_CV_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_CV_adjusted_EAS_boot,R2_adjusted_EAS_boot)

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(beta_RV_validation_raw_EUR,beta_RV_validation_raw_SAS,beta_RV_validation_raw_AMR,beta_RV_validation_raw_AFR,beta_RV_validation_raw_EAS), 
                             beta_se_raw = c(beta_RV_se_validation_raw_EUR,beta_RV_se_validation_raw_SAS,beta_RV_se_validation_raw_AMR,beta_RV_se_validation_raw_AFR,beta_RV_se_validation_raw_EAS), 
                             R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                             R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                             beta_adjusted = c(beta_RV_validation_adjusted_EUR,beta_RV_validation_adjusted_SAS,beta_RV_validation_adjusted_AMR,beta_RV_validation_adjusted_AFR,beta_RV_validation_adjusted_EAS), 
                             beta_se_adjusted = c(beta_RV_se_validation_adjusted_EUR,beta_RV_se_validation_adjusted_SAS,beta_RV_se_validation_adjusted_AMR,beta_RV_se_validation_adjusted_AFR,beta_RV_se_validation_adjusted_EAS), 
                             R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                             R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))

RV_Boot_Results <- data.frame(trait = trait,beta_RV_raw_EUR_boot,R2_raw_EUR_boot,beta_RV_raw_SAS_boot,R2_raw_SAS_boot,
                              beta_RV_raw_AMR_boot,R2_raw_AMR_boot,beta_RV_raw_AFR_boot,R2_raw_AFR_boot,
                              beta_RV_raw_EAS_boot,R2_raw_EAS_boot,beta_RV_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                              beta_RV_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_RV_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                              beta_RV_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_RV_adjusted_EAS_boot,R2_adjusted_EAS_boot)

write.csv(CV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"Best_Betas_RICECV.csv"),row.names = FALSE)
write.csv(CV_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Bootstraps_RICECV.csv"),row.names = FALSE)
write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"Best_Betas_RICERV.csv"),row.names = FALSE)
write.csv(RV_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Bootstraps_RICERV.csv"),row.names = FALSE)