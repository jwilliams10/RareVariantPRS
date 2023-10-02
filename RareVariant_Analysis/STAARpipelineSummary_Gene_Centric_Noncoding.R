##########################################################
# Summarization and visualization of gene-centric 
# noncoding analysis results using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
##########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

###########################################################
#           User Input
###########################################################
## aGDS directory
agds_dir <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/train_agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/nullmodels_staar/Train_Null_Model_LDL.RData"))
## Known loci
known_loci <- NULL

## results path
input_path <- "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/"
output_path <- input_path
## number of jobs
gene_centric_noncoding_jobs_num <- 387
## results name
gene_centric_results_name <- "UKBB_WES_LDL_NonCoding_Train"

## QC_label
QC_label <- "annotation/info/QC_label"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## variant_type
variant_type <- "SNV"
## method_cond
method_cond <- "optimal"
## alpha level
alpha <- 1

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_lipids/Data/agds/train_Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## ncRNA
ncRNA_jobs_num <- 223
ncRNA_input_path <- "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentric_ncRNA/"
ncRNA_output_path <- ncRNA_input_path
ncRNA_results_name <- "UKBB_WES_LDL_ncRNA_Train"

###########################################################
#           Main Function 
###########################################################
## gene info
Gene_Centric_Noncoding_Results_Summary(agds_dir=agds_dir,gene_centric_noncoding_jobs_num=gene_centric_noncoding_jobs_num,
                                       input_path=input_path,output_path=output_path,gene_centric_results_name=gene_centric_results_name,
                                       ncRNA_jobs_num=ncRNA_jobs_num,ncRNA_input_path=ncRNA_input_path,
                                       ncRNA_output_path=ncRNA_output_path,ncRNA_results_name=ncRNA_results_name,
                                       obj_nullmodel=obj_nullmodel,known_loci=known_loci,
                                       method_cond=method_cond,
                                       QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                       alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)