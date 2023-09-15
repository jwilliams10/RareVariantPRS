##########################################################
# Annotate Rare Variants in Coding Masks
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
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(dplyr)

source("/spin1/home/linux/williamsjacr/RareVariantPRS/RareVariant_PRS_Train/Burden_PRS.R")
source("~/RareVariantPRS/RareVariant_PRS_Train/Sliding_Window_Burden_PRS.R")


###########################################################
#           User Input
###########################################################

### Significant Results 
sliding_window_sig <- read_csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/sliding_window_sig.csv")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/train_agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Train_Null_Model_LDL.RData"))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

###########################################################

PRS <- NULL
for(i in 1:nrow(sliding_window_sig)){
  ## Chr
  chr <- sliding_window_sig$Chr[i]
  ## Start
  start_loc <- sliding_window_sig$Start[i]
  ## End
  end_loc <- sliding_window_sig$End[i]
  ## Beta
  BETA <- sliding_window_sig$Burden_Est[i]
  
  ### gds file
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  
  if(i == 1){
    PRS <- Sliding_Window_Burden_PRS(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,
                                     rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                     BETA = BETA,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
  }else{
    PRS[,2] <- PRS[,2] + Sliding_Window_Burden_PRS(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        BETA = BETA,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)[,2]
  }
  seqClose(genofile) 
}

