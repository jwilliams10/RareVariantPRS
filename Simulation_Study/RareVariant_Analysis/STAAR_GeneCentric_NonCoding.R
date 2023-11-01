#####################################################################
# Gene-centric analysis for noncoding rare variants of protein-coding 
# genes using STAARpipeline
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
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(dplyr)

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
output_path <- "/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/"

###########################################################
#           Main Function 
###########################################################
genes_info_chr <- genes_info[genes_info[,2]==22,]

## aGDS file
genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

results_noncoding <- c()
for(kk in 1:nrow(genes_info_chr))
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Noncoding(chr=22,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_noncoding <- append(results_noncoding,results)
}

save(results_noncoding,file=paste0(output_path,"Train_Analysis",i,".Rdata"))

seqClose(genofile)


##################### Summarize

rm(list=setdiff(ls(), "i"))
gc()

load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_Analysis",i,".Rdata"))
noncoding_sig <- NULL
for(j in 1:length(results_noncoding)){
  noncoding_sig <- rbind(noncoding_sig,unlist(results_noncoding[[j]][,c(1:3,47,90)]))
}

noncoding_sig <- as.data.frame(noncoding_sig)
noncoding_sig$Chr <- as.numeric(noncoding_sig$Chr) 

save(noncoding_sig,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/noncoding_sig",i,".Rdata"))









##########################################################
# Annotate Rare Variants in Noncoding Masks
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
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(dplyr)

## source code
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Burden_Effect_Size.R")
source("~/RareVariantPRS/Simulation_Study/RareVariant_PRS/Gene_Centric_Noncoding_Burden_Effect_Size_Jake.R")

###########################################################
#           User Input
###########################################################

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












rm(list=setdiff(ls(), "i"))

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricNoncoding/chr22.Rdata")
noncoding_sig <- NULL
for(j in 1:length(results_noncoding)){
  noncoding_sig <- rbind(noncoding_sig,unlist(results_noncoding[[j]][,1:3]))
}

noncoding_sig <- as.data.frame(noncoding_sig)
noncoding_sig$Chr <- as.numeric(noncoding_sig$Chr)
noncoding_sig_overall <- noncoding_sig
colnames(noncoding_sig_overall)[1] <- "Gene"

noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes",i,".csv"))
noncoding_sig <- left_join(noncoding_sig_overall,noncoding_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Tune_Null_Model",i,".RData")))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Noncoding.RData")

G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[,!is.na(noncoding_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- noncoding_sig$Burden_pvalue
  }else{
    p_values <- noncoding_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_gene_centric_noncoding))
    beta[p_values < thresholds[k]] <- noncoding_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_gene_centric_noncoding%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/tune_prs_mat",i,".RData"))













rm(list=setdiff(ls(), "i"))

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricNoncoding/chr22.Rdata")
noncoding_sig <- NULL
for(j in 1:length(results_noncoding)){
  noncoding_sig <- rbind(noncoding_sig,unlist(results_noncoding[[j]][,1:3]))
}

noncoding_sig <- as.data.frame(noncoding_sig)
noncoding_sig$Chr <- as.numeric(noncoding_sig$Chr)
noncoding_sig_overall <- noncoding_sig
colnames(noncoding_sig_overall)[1] <- "Gene"

noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes",i,".csv"))
noncoding_sig <- left_join(noncoding_sig_overall,noncoding_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Validation_Null_Model",i,".RData")))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Noncoding.RData")

G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[,!is.na(noncoding_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- noncoding_sig$Burden_pvalue
  }else{
    p_values <- noncoding_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_gene_centric_noncoding))
    beta[p_values < thresholds[k]] <- noncoding_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_gene_centric_noncoding%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/validation_prs_mat",i,".RData"))




