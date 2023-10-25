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
library(readr)
library(dplyr)

## source code
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Burden_Effect_Size.R")
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Sliding_Window_Burden_Effect_Size.R")

###########################################################
#           User Input
###########################################################

i <- as.numeric(commandArgs(TRUE)[1])

### Significant Results 
load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/slidingwindow_sig",i,".Rdata"))
colnames(slidingwindow_sig) <- c("Chr","Start","End","Burden_1_1","STAAR_O")

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
for(j in 1:nrow(slidingwindow_sig)){
  ## Chr
  chr <- slidingwindow_sig$Chr[j]
  ## Start
  start_loc <- slidingwindow_sig$Start[j]
  ## End
  end_loc <- slidingwindow_sig$End[j]
  
  a <- Sliding_Window_Burden_Effect_Size(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,
                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
  a <- data.frame(Chr = chr,Start = start_loc,End = end_loc,Number_SNV = a[[4]],cMAC = a[[5]],Burden_Score_Stat = a[[6]],Burden_SE_Score = a[[7]],Burden_pvalue = a[[8]],Burden_Est = a[[9]], Burden_SE_Est = a[[10]])
  effect_sizes <- rbind(effect_sizes,a)
}
seqClose(genofile) 

slidingwindow_sig <- inner_join(slidingwindow_sig,effect_sizes)

write.csv(slidingwindow_sig,row.names = FALSE,file = paste0(output_path <- "/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes",i,".csv"))


		
