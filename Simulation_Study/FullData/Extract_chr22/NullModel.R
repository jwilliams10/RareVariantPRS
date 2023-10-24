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

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sample_phenotype.RData")
colnames(phenotype) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

obj.STAAR.UKB.LDL <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData")