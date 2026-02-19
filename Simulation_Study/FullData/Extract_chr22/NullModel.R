rm(list = ls())
# =============================================================================
# [RICE-ANNOTATION] FullData/Extract_chr22/NullModel.R
# Purpose: Fits a STAAR null model on a synthetic phenotype for the extracted chr22 dataset (used as a pre-step for association scans).
#
# Paper linkage:
#   - Manuscript; Methods → Simulation Study; Results → Simulation Study Results (Fig. 3).
#   - Supplementary Data: Supplementary Data 1 (sheet “S1 Sample Sizes; Sim. Study”).
#   - Supplementary Figures: Supp. Fig. 1–2; Supp. Note → “Ancestry Adjusted PRS”.
# =============================================================================

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
