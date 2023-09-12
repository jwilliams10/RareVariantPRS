##########################################################
# Annotate Rare Variants in Sliding Windows
# using STAARpipelineSummary
# Xihao Li, Zilin Li
# 06/29/2023
##########################################################

rm(list=ls())
gc()

## load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # or library(Matrix)

## source code
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/Burden_Effect_Size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/Sliding_Window_Burden_Effect_Size_silent.R")

###########################################################
#           User Input
###########################################################
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/obj.STAAR.CAD.freeze9.noage.20230620.Rdata"))

## Parameter
QC_label <- "annotation/filter"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

###########################################################
#      1-41003016-41005015
###########################################################

## Chr
chr <- 1
## location
start_loc <- 41003016
end_loc <- 41005015

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Sliding_Window_Burden_Effect_Size(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
								
seqClose(genofile)




		
