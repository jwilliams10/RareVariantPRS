# Step3: Run NULL model
rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)


## Phenotype
phenotype <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/UKB_WES_lipids/tmp.TC.20211014.Rdata"))
objects()

## test
lm.test <- lm(TCadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype)
summary(lm.test)

## load GRM
sgrm <- get(load("/n/holystore01/LABS/xlin/Lab/rdey/UKB_Analysis/FastSparseGRM/grmout/grmout_nRandomSNPs_0.sGRM.RData"))
sample_id <- unlist(lapply(strsplit(colnames(sgrm),"_"),`[[`,1))

colnames(sgrm) <- sample_id
rownames(sgrm) <- sample_id


a <- Sys.time()
### fit null model
obj.STAAR.UKB.TC <- fit_null_glmmkin(TCadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype, kins = sgrm, kins_cutoff = 0.022, id = "userId", use_sparse = TRUE,family = gaussian(link = "identity"), verbose=T)
b <- Sys.time()
b - a

objects()

save(obj.STAAR.UKB.TC,file = "/n/holystore01/LABS/xlin/Lab/xihao_zilin/UKB_WES_lipids/obj.STAAR.UKB.TC.20211014.Rdata")







