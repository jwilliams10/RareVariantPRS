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
colnames(phenotype) <- c("IID","FID","Y")

obj.STAAR.UKB.LDL <- fit_nullmodel(Y~1, data = phenotype,id = "IID",kins = NULL,family = gaussian(link = "identity"))

save(obj.STAAR.UKB.LDL,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData")
