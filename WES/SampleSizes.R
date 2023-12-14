rm(list = ls())

pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

rm(pheno_vad)

samps <- NULL

for(trait in c("Asthma","CAD","T2D","Breast","Prostate","BMI","TC","LDL","HDL","logTG","Height")){
  y <- pheno_train[,trait]
  
  n_train <- sum(!is.na(y))
  n_train_control <- sum(y == 0,na.rm = TRUE)
  n_train_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_tuning[,trait]
  
  n_tune <- sum(!is.na(y))
  n_tune_control <- sum(y == 0,na.rm = TRUE)
  n_tune_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_EUR[,trait]
  
  n_vad_EUR <- sum(!is.na(y))
  n_vad_EUR_control <- sum(y == 0,na.rm = TRUE)
  n_vad_EUR_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_NonEur[,trait]
  
  n_vad_NonEur <- sum(!is.na(y))
  n_vad_NonEur_control <- sum(y == 0,na.rm = TRUE)
  n_vad_NonEur_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_UNK[,trait]
  
  n_vad_UNK <- sum(!is.na(y))
  n_vad_UNK_control <- sum(y == 0,na.rm = TRUE)
  n_vad_UNK_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_EAS[,trait]
  
  n_vad_EAS <- sum(!is.na(y))
  n_vad_EAS_control <- sum(y == 0,na.rm = TRUE)
  n_vad_EAS_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_SAS[,trait]
  
  n_vad_SAS <- sum(!is.na(y))
  n_vad_SAS_control <- sum(y == 0,na.rm = TRUE)
  n_vad_SAS_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_MIX[,trait]
  
  n_vad_MIX <- sum(!is.na(y))
  n_vad_MIX_control <- sum(y == 0,na.rm = TRUE)
  n_vad_MIX_cases <- sum(y == 1,na.rm = TRUE)
  
  y <- pheno_vad_AFR[,trait]
  
  n_vad_AFR <- sum(!is.na(y))
  n_vad_AFR_control <- sum(y == 0,na.rm = TRUE)
  n_vad_AFR_cases <- sum(y == 1,na.rm = TRUE)
  
  samps_temp <- data.frame(Trait = rep(trait,7),Ancestry = c("EUR","NONEUR","UNK","EAS","SAS","MIX","AFR"), Train_Total = c(n_train,0,0,0,0,0,0), Train_Controls = c(n_train_control,0,0,0,0,0,0), Train_Cases = c(n_train_cases,0,0,0,0,0,0),
                           Tune_Total = c(n_tune,0,0,0,0,0,0), Tune_Controls = c(n_tune_control,0,0,0,0,0,0), Tune_Cases = c(n_tune_cases,0,0,0,0,0,0),
                           Validation_Total = c(n_vad_EUR,n_vad_NonEur,n_vad_UNK,n_vad_EAS,n_vad_SAS,n_vad_MIX,n_vad_AFR), 
                           Validation_Controls = c(n_vad_EUR_control,n_vad_NonEur_control,n_vad_UNK_control,n_vad_EAS_control,n_vad_SAS_control,n_vad_MIX_control,n_vad_AFR_control), 
                           Validation_Cases = c(n_vad_EUR_cases,n_vad_NonEur_cases,n_vad_UNK_cases,n_vad_EAS_cases,n_vad_SAS_cases,n_vad_MIX_cases,n_vad_AFR_cases))
  
  samps <- rbind(samps,samps_temp)
  
}

rm(list=setdiff(ls(), "samps"))

samps$Train_Controls[samps$Trait %in% c("BMI","TC","LDL","HDL","logTG","Height")] <- 0
samps$Train_Cases[samps$Trait %in% c("BMI","TC","LDL","HDL","logTG","Height")] <- 0
samps$Tune_Controls[samps$Trait %in% c("BMI","TC","LDL","HDL","logTG","Height")] <- 0
samps$Tune_Cases[samps$Trait %in% c("BMI","TC","LDL","HDL","logTG","Height")] <- 0
samps$Validation_Controls[samps$Trait %in% c("BMI","TC","LDL","HDL","logTG","Height")] <- 0
samps$Validation_Cases[samps$Trait %in% c("BMI","TC","LDL","HDL","logTG","Height")] <- 0
