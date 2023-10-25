##########################################################
# Annotate Rare Variants in Noncoding Masks
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

## source code
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Burden_Effect_Size.R")
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Gene_Centric_Noncoding_Burden_Effect_Size_Jake.R")

###########################################################
#           User Input
###########################################################

i <- as.numeric(commandArgs(TRUE)[1])

### Significant Results 
load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/noncoding_sig",i,".Rdata"))
colnames(noncoding_sig) <- c("Gene","Chr","Category","Burden_1_1","STAAR_O")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Train_Null_Model",i,".RData")))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_Annotation_name_catalog.Rdata"))

genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

effect_sizes <- NULL
for(j in 1:nrow(noncoding_sig)){
  ## Chr
  chr <- noncoding_sig$Chr[j]
  ## Gene name
  gene_name <- noncoding_sig$Gene[j]
  ## Coding mask
  category <- noncoding_sig$Category[j]
  
  a <- Gene_Centric_Noncoding_Burden_Effect_Size_Jake(chr=chr,gene_name=gene_name,category=category ,
                                                 genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
  a <- data.frame(Gene = gene_name,Chr = chr,Category = category,Number_SNV = a[[4]],cMAC = a[[5]],Burden_Score_Stat = a[[6]],Burden_SE_Score = a[[7]],Burden_pvalue = a[[8]],Burden_Est = a[[9]], Burden_SE_Est = a[[10]])
  effect_sizes <- rbind(effect_sizes,a)
}

seqClose(genofile) 

noncoding_sig <- inner_join(noncoding_sig,effect_sizes)

write.csv(noncoding_sig,row.names = FALSE,file = paste0(output_path <- "/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes",i,".csv"))
