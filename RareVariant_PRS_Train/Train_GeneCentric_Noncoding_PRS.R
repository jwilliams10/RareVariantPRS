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
source("~/RareVariantPRS/RareVariant_PRS_Train/Gene_Centric_Noncoding_Burden_PRS.R")


###########################################################
#           User Input
###########################################################

### Significant Results 
noncoding_sig <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/noncoding_sig.csv")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/train_agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Train_Null_Model_LDL.RData"))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/train_Annotation_name_catalog.Rdata"))

PRS <- NULL
for(i in 1:nrow(noncoding_sig)){
  ## Chr
  chr <- noncoding_sig$Chr[i]
  ## Gene name
  gene_name <- noncoding_sig$Gene[i]
  ## Coding mask
  category <- noncoding_sig$Category[i]
  ## Beta
  BETA <- noncoding_sig$Burden_Est[i]
  
  ### gds file
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  
  if(i == 1){
    PRS <- Gene_Centric_Noncoding_Burden_PRS(chr=chr,gene_name=gene_name,category=category ,
                                          genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                          BETA = BETA,
                                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
  }else{
    PRS[,2] <- PRS[,2] + Gene_Centric_Noncoding_Burden_PRS(chr=chr,gene_name=gene_name,category=category ,
                                                        genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        BETA = BETA,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)[,2]
  }
  seqClose(genofile) 
}

