##########################################################
# Manually create a QC_label with all "PASS" (for Xinan)
# Xihao Li, Zilin Li, Hufeng Zhou
# 10/30/2021
##########################################################

rm(list=ls())
gc()

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

QC_label <- list()

chr <- 22

print(chr)
### gds file
gds.path <- "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds"

genofile <- seqOpen(gds.path, readonly = FALSE)

F_Missing <- seqGetData(genofile, "annotation/info/F_MISSING")
QC_label[[chr]] <- ifelse(F_Missing<0.1,"PASS","FAIL") 

sum(QC_label[[chr]]=="FAIL")

Anno.folder <- index.gdsn(genofile, "annotation/info")
add.gdsn(Anno.folder, "QC_label", val=QC_label[[chr]], compress="LZMA_ra", closezip=TRUE)

seqClose(genofile)







