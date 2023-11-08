#####################################################################
# Gene-centric analysis for coding rare variants using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
#####################################################################

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

###########################################################
#           User Input
###########################################################

i <- as.numeric(commandArgs(TRUE)[1])

## job nums
jobs_num <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Train_Null_Model",i,".RData")))

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

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
output_path <- "/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/"

###########################################################
#           Main Function 
###########################################################

genes_info_chr <- genes_info[genes_info[,2]==22,]


### gds file
genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

results_coding <- c()

for(kk in 1:nrow(genes_info_chr)){
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Coding(chr=22,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_coding <- append(results_coding,results)
}

save(results_coding,file=paste0(output_path,"Train_Analysis",i,".Rdata"))

seqClose(genofile)


##################### Summarize

rm(list=setdiff(ls(), "i"))
gc()

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/Train_Analysis",i,".Rdata"))
coding_sig <- NULL
for(j in 1:length(results_coding)){
  coding_sig <- rbind(coding_sig,unlist(results_coding[[j]][,c(1:3,47,90)]))
}

coding_sig <- as.data.frame(coding_sig)
coding_sig$Chr <- as.numeric(coding_sig$Chr) 

save(coding_sig,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/coding_sig",i,".Rdata"))






rm(list=setdiff(ls(), "i"))
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
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Gene_Centric_Coding_Burden_Effect_Size_Jake.R")

###########################################################
#           User Input
###########################################################

### Significant Results 
load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/coding_sig",i,".Rdata"))
colnames(coding_sig) <- c("Gene","Chr","Category","Burden_1_1","STAAR_O")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Train_Null_Model",i,".RData")))

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
# effect_sizes_Jake <- NULL
for(j in 1:nrow(coding_sig)){
  ## Chr
  chr <- coding_sig$Chr[j]
  ## Gene name
  gene_name <- coding_sig$Gene[j]
  ## Coding mask
  category <- coding_sig$Category[j]
  
  a <- Gene_Centric_Coding_Burden_Effect_Size_Jake(chr=chr,gene_name=gene_name,category=category ,
                                                   genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                   QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
  a <- data.frame(Gene = gene_name,Chr = chr,Category = category,Number_SNV = a[[4]],cMAC = a[[5]],Burden_Score_Stat = a[[6]],Burden_SE_Score = a[[7]],Burden_pvalue = a[[8]],Burden_Est = a[[9]], Burden_SE_Est = a[[10]])
  effect_sizes <- rbind(effect_sizes,a)
}
seqClose(genofile) 

coding_sig <- inner_join(coding_sig,effect_sizes)

write.csv(coding_sig,row.names = FALSE,file = paste0(output_path <- "/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/Train_EffectSizes",i,".csv"))




rm(list=setdiff(ls(), "i"))

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricCoding/chr22.Rdata")
coding_sig <- NULL
for(j in 1:length(results_coding)){
  coding_sig <- rbind(coding_sig,unlist(results_coding[[j]][,1:3]))
}

coding_sig <- as.data.frame(coding_sig)
coding_sig$Chr <- as.numeric(coding_sig$Chr)
coding_sig_overall <- coding_sig
colnames(coding_sig_overall)[1] <- "Gene"

coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/Train_EffectSizes",i,".csv"))
coding_sig <- left_join(coding_sig_overall,coding_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Tune_Null_Model",i,".RData")))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Coding.RData")

G_star_gene_centric_coding <- G_star_gene_centric_coding[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_gene_centric_coding <- G_star_gene_centric_coding[,!is.na(coding_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- coding_sig$Burden_pvalue
  }else{
    p_values <- coding_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_gene_centric_coding))
    beta[p_values < thresholds[k]] <- coding_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_gene_centric_coding%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/tune_prs_mat",i,".RData"))



rm(list=setdiff(ls(), "i"))

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricCoding/chr22.Rdata")
coding_sig <- NULL
for(j in 1:length(results_coding)){
  coding_sig <- rbind(coding_sig,unlist(results_coding[[j]][,1:3]))
}

coding_sig <- as.data.frame(coding_sig)
coding_sig$Chr <- as.numeric(coding_sig$Chr)
coding_sig_overall <- coding_sig
colnames(coding_sig_overall)[1] <- "Gene"

coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/Train_EffectSizes",i,".csv"))
coding_sig <- left_join(coding_sig_overall,coding_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/nullmodels_staar/Validation_Null_Model",i,".RData")))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Coding.RData")

G_star_gene_centric_coding <- G_star_gene_centric_coding[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_gene_centric_coding <- G_star_gene_centric_coding[,!is.na(coding_sig$Burden_1_1)]

### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)
coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- coding_sig$Burden_pvalue
  }else{
    p_values <- coding_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_gene_centric_coding))
    beta[p_values < thresholds[k]] <- coding_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_gene_centric_coding%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/validation_prs_mat",i,".RData"))


