rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)
library(SCANG)
library(STAARpipeline)

i <- as.numeric(commandArgs(TRUE)[1])

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Train_All",i,".txt"))

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Train_All",i,".txt")))

pheno_train <- inner_join(pheno_train,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_nullmodel(Y~PRS, data = pheno_train,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Train_Null_Model",i,".RData"))


rm(list=setdiff(ls(), "i"))

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")
pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Tune_All",i,".txt"))

pheno_tune <- inner_join(pheno_tune,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_nullmodel(Y~PRS, data = pheno_tune,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Tune_Null_Model",i,".RData"))


rm(list=setdiff(ls(), "i"))

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")

common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Validation_All",i,".txt"))

pheno_vad <- inner_join(pheno_vad,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_nullmodel(Y~PRS, data = pheno_vad,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Validation_Null_Model",i,".RData"))
  
