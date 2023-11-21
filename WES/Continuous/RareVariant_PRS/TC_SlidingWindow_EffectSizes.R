##########################################################
# Annotate Rare Variants in Sliding Windows
# using STAARpipelineSummary
# Xihao Li, Zilin Li
# 06/29/2023
##########################################################

rm(list=ls())
gc()

trait <- "TC"

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
source("~/RareVariantPRS/WES/Continuous/RareVariant_PRS/Burden_Effect_Size.R")
source("~/RareVariantPRS/WES/Continuous/RareVariant_PRS/Sliding_Window_Burden_Effect_Size.R")

###########################################################
#           User Input
###########################################################
### Significant Results 
sliding_window_sig <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_sliding_window_sig.csv"))
sliding_window_sig <- sliding_window_sig[,c(1,2,3,4,5,6,48,62,91)]
colnames(sliding_window_sig) <- c("IDK","Chr","Start","End","Number_SNV","SKAT_1_25","Burden_1_1","ACAT_V_1_25","STAAR_O")

## aGDS directory
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Train_Null_Model.RData")))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

###########################################################

arrayid <- as.numeric(commandArgs(TRUE)[1])

sliding_window_sig <- sliding_window_sig[sliding_window_sig$Chr == arrayid,]

effect_sizes <- NULL
for(i in 1:nrow(sliding_window_sig)){
  ## Chr
  chr <- sliding_window_sig$Chr[i]
  ## Start
  start_loc <- sliding_window_sig$Start[i]
  ## End
  end_loc <- sliding_window_sig$End[i]
  
  ### gds file
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  
  a <- Sliding_Window_Burden_Effect_Size(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,
                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
  a <- data.frame(IDK = row.names(a),Chr = chr,Start = start_loc,End = end_loc,Number_SNV = a[[4]],cMAC = a[[5]],Burden_Score_Stat = a[[6]],Burden_SE_Score = a[[7]],Burden_pvalue = a[[8]],Burden_Est = a[[9]], Burden_SE_Est = a[[10]])
  effect_sizes <- rbind(effect_sizes,a)
  seqClose(genofile) 
}

sliding_window_sig <- sliding_window_sig[,-c(1)]
effect_sizes <- effect_sizes[,-c(1)]

sliding_window_sig <- inner_join(sliding_window_sig,effect_sizes)

write.csv(sliding_window_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_sliding_window_sig_chr",arrayid,".csv"))


		
