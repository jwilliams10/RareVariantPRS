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

pheno_train <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt")
colnames(pheno_train) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Train_All.txt")

pheno_train <- inner_join(pheno_train,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + prs, data = pheno_train,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = "/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Train_Null_Model_LDL.RData")

# obj_nullmodel_SCANG_STAAR <- staar2scang_nullmodel(obj.STAAR.UKB.LDL)
# 
# save(obj_nullmodel_SCANG_STAAR,file="/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Train_Null_Model_SCANG_LDL.RData")

rm(list = ls())

pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
colnames(pheno_tune) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Tune_All.txt")

pheno_tune <- inner_join(pheno_tune,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + prs, data = pheno_tune,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = "/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Tune_Null_Model_LDL.RData")

# obj_nullmodel_SCANG_STAAR <- staar2scang_nullmodel(obj.STAAR.UKB.LDL)
# 
# save(obj_nullmodel_SCANG_STAAR,file="/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Tune_Null_Model_SCANG_LDL.RData")

rm(list = ls())

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt")

pheno_vad <- inner_join(pheno_vad,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + prs, data = pheno_vad,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = "/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Validation_Null_Model_LDL.RData")

# obj_nullmodel_SCANG_STAAR <- staar2scang_nullmodel(obj.STAAR.UKB.LDL)
# 
# save(obj_nullmodel_SCANG_STAAR,file="/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Validation_Null_Model_SCANG_LDL.RData")
