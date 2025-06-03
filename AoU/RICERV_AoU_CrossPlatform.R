rm(list = ls())
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(dplyr)
library(RISCA)
library(boot)
library(stringr)
library(caret)
library(ranger)
library(glmnet)

source("~/RareVariantPRS/Simulation_Study/FullData/G_star/g_star_gene_centric_coding.R")

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")

trait <- c(continuous_traits)[as.numeric(commandArgs(TRUE)[1])]

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
sampleids_all <- read.table("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega.fam", quote="\"", comment.char="")
unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")
ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids_all[,2],]
ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% unrels_nRandomSNPs_0[,2],]

final_coef <- read.csv(paste0("/data/williamsjacr/AoU_Results/",trait,"_final_coef.csv"))
final_coef <- final_coef[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","Beta")]

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))

## Null Model
obj_nullmodel_train <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Train_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(ukb_pheno$IID)

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/Annotation_name_catalog.Rdata"))

G_star_gene_centric_coding <- list()

for(i in 1:nrow(final_coef)){
  ## Chr
  chr <- final_coef$Chr[i]
  ## Gene name
  gene_name <- final_coef$Gene[i]
  ## Coding mask
  category <- final_coef$Category[i]
  
  ### gds file
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  
  G_star_gene_centric_coding[[i]] <- tryCatch(Gene_Centric_Coding_G_Star(chr=chr,gene_name=gene_name,category=category ,
                                                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE),error = function(cond) {return(matrix(0,ncol = 1,nrow = length(obj_nullmodel$id_include)))})
  seqClose(genofile) 
  G_star_gene_centric_coding[[i]] <- final_coef$Beta[i]*G_star_gene_centric_coding[[i]]
} 

G_star_gene_centric_coding <- do.call(cbind,G_star_gene_centric_coding)

RICE_RV_PRS <- data.frame(IID = obj_nullmodel$id_include,PRS = rowSums(G_star_gene_centric_coding))
colnames(RICE_RV_PRS) <- c("IID","PRS")

write.csv(RICE_RV_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_RICERV_PRS.csv"),row.names = FALSE)
