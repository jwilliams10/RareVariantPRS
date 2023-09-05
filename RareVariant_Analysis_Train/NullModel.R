rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)

pheno_train <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt")
colnames(pheno_train) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

common_prs <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Train.txt")

pheno_train <- inner_join(pheno_train,common_prs,by = "IID")

obj.STAAR.UKB.LDL <- fit_null_glm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + prs, data = pheno_train,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = "/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Train_Null_Model_LDL.RData")