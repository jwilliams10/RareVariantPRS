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

source("/spin1/home/linux/williamsjacr/RareVariantPRS/RareVariant_PRS/Burden_PRS.R")
source("~/RareVariantPRS/RareVariant_PRS/Gene_Centric_Noncoding_Burden_PRS.R")


###########################################################
#           User Input
###########################################################

### Significant Results 
Train_Effect_Sizes_All <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Train_Effect_Sizes_All.csv")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/tune_agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Tune_Null_Model_LDL.RData"))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/tune_Annotation_name_catalog.Rdata"))
thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

PRS <- NULL

arrayid <- as.numeric(commandArgs(TRUE)[1])


if(arrayid>330){
  Burden <- 1
  arrayid <- arrayid - 330
}else{
  Burden <- 0
}

if(arrayid <= 22){
  threshold <- 1
}else if(arrayid <= 44){
  threshold <- 2
}else if(arrayid <= 66){
  threshold <- 3
}else if(arrayid <= 88){
  threshold <- 4
}else if(arrayid <= 110){
  threshold <- 5
}else if(arrayid <= 132){
  threshold <- 6
}else if(arrayid <= 154){
  threshold <- 7
}else if(arrayid <= 176){
  threshold <- 8
}else if(arrayid <= 198){
  threshold <- 9
}else if(arrayid <= 220){
  threshold <- 10
}else if(arrayid <= 242){
  threshold <- 11
}else if(arrayid <= 264){
  threshold <- 12
}else if(arrayid <= 286){
  threshold <- 13
}else if(arrayid <= 308){
  threshold <- 14
}else{
  threshold <- 15
}

chr <- arrayid %% 22
if(chr == 0){chr <- 22}

Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$Chr == chr,]

if(Burden == 0){
  Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$STAAR_O <= thresholds[threshold],]
}else{
  Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$Burden_1_1 <= thresholds[threshold],]
}

if(nrow(Train_Effect_Sizes_All) == 0){
  PRS <- data.frame(ID = 1:length(obj_nullmodel$id_include),PRS = 0)
}else{
  for(i in 1:nrow(Train_Effect_Sizes_All)){
    ## Chr
    chr <- Train_Effect_Sizes_All$Chr[i]
    ## Gene name
    gene_name <- Train_Effect_Sizes_All$Gene[i]
    ## Coding mask
    category <- Train_Effect_Sizes_All$Category[i]
    ## Beta
    BETA <- Train_Effect_Sizes_All$Burden_Est[i]
    
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
}

write.csv(PRS,file = paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Tune_PRS_Array_",arrayid,".csv"),row.names = FALSE)
