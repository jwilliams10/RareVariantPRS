#####################################################################
# Sliding window analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
#####################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

i <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           User Input
###########################################################
## job nums
jobs_num <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Train_Null_Model",i,".RData")))

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## sliding_window_length
sliding_window_length <- 2000

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/"

###########################################################
#           Main Function 
###########################################################

genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

start_loc <- jobs_num$start_loc[1]
end_loc <- start_loc + (sliding_window_length/2)*20 - 1
end_loc_final <- jobs_num$end_loc[1]

results_sliding_window <- c()
for(kk in 1:(((end_loc_final - start_loc)/((sliding_window_length/2)*20)) + 1)){
  print(kk)
  start_loc_sub <- start_loc + (sliding_window_length/2)*20*(kk-1)
  end_loc_sub <- end_loc + (sliding_window_length/2)*20*(kk-1) + (sliding_window_length/2)
  
  end_loc_sub <- min(end_loc_sub,jobs_num$end_loc[1])
  
  results <- c()
  if(start_loc_sub < end_loc_sub){
    results <- try(Sliding_Window(chr=22,start_loc=start_loc_sub,end_loc=end_loc_sub,
                                  sliding_window_length=sliding_window_length,type="multiple",
                                  genofile=genofile,obj_nullmodel=obj_nullmodel,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
    
    if(class(results)[1]!="try-error")
    {
      results_sliding_window <- rbind(results_sliding_window,results)
    }
    
  }
}

save(results_sliding_window,file=paste0(output_path,"Train_Analysis",i,".Rdata"))

seqClose(genofile)

##################### Summarize

rm(list=setdiff(ls(), "i"))
gc()

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_Analysis",i,".Rdata"))

slidingwindow_sig <- data.frame(Chr = unlist(results_sliding_window[,1]),Start = unlist(results_sliding_window[,2]),End = unlist(results_sliding_window[,3]),Burden_1_1 = unlist(results_sliding_window[,47]),STAARO = unlist(results_sliding_window[,90]))

slidingwindow_sig$Chr <- as.numeric(slidingwindow_sig$Chr) 

save(slidingwindow_sig,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/slidingwindow_sig",i,".Rdata"))






##########################################################
# Annotate Rare Variants in Sliding Windows
# using STAARpipelineSummary
# Xihao Li, Zilin Li
# 06/29/2023
##########################################################

rm(list=setdiff(ls(), "i"))
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





rm(list=setdiff(ls(), "i"))

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/SlidingWindow/chr22.Rdata")
slidingwindow_sig <- data.frame(Chr = unlist(results_sliding_window[,1]),Start = unlist(results_sliding_window[,2]),End = unlist(results_sliding_window[,3]))
slidingwindow_sig$Chr <- as.numeric(slidingwindow_sig$Chr)
slidingwindow_sig_overall <- slidingwindow_sig

slidingwindow_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes",i,".csv"))
slidingwindow_sig <- left_join(slidingwindow_sig_overall,slidingwindow_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Tune_Null_Model",i,".RData")))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/SlidingWindow.RData")

G_star_sliding_window <- G_star_sliding_window[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_sliding_window <- G_star_sliding_window[,!is.na(slidingwindow_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

slidingwindow_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- slidingwindow_sig$Burden_pvalue
  }else{
    p_values <- slidingwindow_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_sliding_window))
    beta[p_values < thresholds[k]] <- slidingwindow_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_sliding_window%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/tune_prs_mat",i,".RData"))







rm(list=setdiff(ls(), "i"))

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/SlidingWindow/chr22.Rdata")
slidingwindow_sig <- data.frame(Chr = unlist(results_sliding_window[,1]),Start = unlist(results_sliding_window[,2]),End = unlist(results_sliding_window[,3]))
slidingwindow_sig$Chr <- as.numeric(slidingwindow_sig$Chr)
slidingwindow_sig_overall <- slidingwindow_sig

slidingwindow_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes",i,".csv"))
slidingwindow_sig <- left_join(slidingwindow_sig_overall,slidingwindow_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Validation_Null_Model",i,".RData")))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/SlidingWindow.RData")

G_star_sliding_window <- G_star_sliding_window[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_sliding_window <- G_star_sliding_window[,!is.na(slidingwindow_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

slidingwindow_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- slidingwindow_sig$Burden_pvalue
  }else{
    p_values <- slidingwindow_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_sliding_window))
    beta[p_values < thresholds[k]] <- slidingwindow_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_sliding_window%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/validation_prs_mat",i,".RData"))