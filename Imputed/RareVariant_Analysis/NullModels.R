rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)
library(SCANG)
library(STAARpipeline)

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")

binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")

for(trait in c(continuous_traits,binary_traits)){
  time <- system.time({
    if(trait %in% continuous_traits){
      
      pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
      common_prs <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Train.csv"))
      
      pheno_train <- inner_join(pheno_train,common_prs,by = "IID")
      
      obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_train,id = "IID",kins = NULL,family = gaussian(link = "identity"))
      save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Train_Null_Model.RData"))
      
      pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
      common_prs <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"))
      
      pheno_tune <- inner_join(pheno_tune,common_prs,by = "IID")
      
      obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_tune,id = "IID",kins = NULL,family = gaussian(link = "identity"))
      save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Tune_Null_Model.RData"))
      
      pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
      common_prs <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
      
      pheno_validation <- inner_join(pheno_validation,common_prs,by = "IID")
      
      obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_validation,id = "IID",kins = NULL,family = gaussian(link = "identity"))
      save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Validation_Null_Model.RData"))
      
    }else{
      pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
      common_prs <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Train.csv"))
      
      pheno_train <- inner_join(pheno_train,common_prs,by = "IID")
      if(trait %in% c("Breast","Prostate")){
        obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_train,id = "IID",kins = NULL,family = binomial(link = "logit"),use_SPA = TRUE)
      }else{
        obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_train,id = "IID",kins = NULL,family = binomial(link = "logit"),use_SPA = TRUE) 
      }
      save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Train_Null_Model.RData"))
      
      pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
      common_prs <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"))
      
      pheno_tune <- inner_join(pheno_tune,common_prs,by = "IID")
      if(trait %in% c("Breast","Prostate")){
        obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_tune,id = "IID",kins = NULL,family = binomial(link = "logit"),use_SPA = TRUE)
      }else{
        obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_tune,id = "IID",kins = NULL,family = binomial(link = "logit"),use_SPA = TRUE) 
      }
      save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Tune_Null_Model.RData"))
      
      pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
      common_prs <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
      
      pheno_validation <- inner_join(pheno_validation,common_prs,by = "IID")
      if(trait %in% c("Breast","Prostate")){
        obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_validation,id = "IID",kins = NULL,family = binomial(link = "logit"),use_SPA = TRUE)
      }else{
        obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + PRS")), data = pheno_validation,id = "IID",kins = NULL,family = binomial(link = "logit"),use_SPA = TRUE) 
      }
      save(obj.STAAR.UKB,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Validation_Null_Model.RData"))
    }
    
  })[3]
  save(time,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/NullModel_Time_",trait,".RData"))
}