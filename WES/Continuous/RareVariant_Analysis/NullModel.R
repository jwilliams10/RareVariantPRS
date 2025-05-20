rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)
library(SCANG)
library(STAARpipeline)

trait <- "BMI"

for(trait in c("BMI","TC","LDL","HDL","logTG","Height")){
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Train_All.txt"))
  
  pheno_train <- inner_join(pheno_train,common_prs,by = "IID")
  
  obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + prs")), data = pheno_train,id = "IID",kins = NULL,family = gaussian(link = "identity"))
  
  save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Train_Null_Model.RData"))
  
  
  
  
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"))
  
  pheno_tune <- inner_join(pheno_tune,common_prs,by = "IID")
  
  obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + prs")), data = pheno_tune,id = "IID",kins = NULL,family = gaussian(link = "identity"))
  
  save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Tune_Null_Model.RData"))
  
  
  
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  
  common_prs <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
  
  pheno_validation <- inner_join(pheno_validation,common_prs,by = "IID")
  
  obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + prs")), data = pheno_validation,id = "IID",kins = NULL,family = gaussian(link = "identity"))
  
  save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Validation_Null_Model.RData"))
  
  
}