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

chr <- as.numeric(commandArgs(TRUE)[1])

print(chr)
### gds file
gds.dir <- "/data/williamsjacr/UKB_WES_Full_Processed_Data/gds/"
gds.path <- paste0(gds.dir,"chr",chr,".gds")

genofile <- seqOpen(gds.path, readonly = FALSE)

F_Missing <- seqGetData(genofile, "annotation/info/F_MISSING")
QC_label <- ifelse(F_Missing<0.1,"PASS","FAIL") 

sum(QC_label=="FAIL")

Anno.folder <- index.gdsn(genofile, "annotation/info")
add.gdsn(Anno.folder, "QC_label", val=QC_label, compress="LZMA_ra", closezip=TRUE)

seqClose(genofile)

