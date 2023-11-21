##########################################################
# Annotate Rare Variants in Noncoding Masks
# using STAARpipelineSummary
# Xihao Li, Zilin Li
# 06/29/2023
##########################################################

rm(list=ls())
gc()

trait <- "LDL"

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
source("~/RareVariantPRS/WES/Continuous/RareVariant_PRS/Burden_Effect_Size.R")
source("~/RareVariantPRS/WES/Continuous//RareVariant_PRS/Gene_Centric_Noncoding_Burden_Effect_Size_Jake.R")

###########################################################
#           User Input
###########################################################
### Significant Results 
noncoding_sig <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_noncoding_sig.csv"))
noncoding_sig <- noncoding_sig[,c(1,2,3,4,5,6,48,62,91)]
colnames(noncoding_sig) <- c("IDK","Gene","Chr","Category","Number_SNV","SKAT_1_25","Burden_1_1","ACAT_V_1_25","STAAR_O")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Train_Null_Model.RData")))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/Annotation_name_catalog.Rdata"))

###########################################################

arrayid <- as.numeric(commandArgs(TRUE)[1])

noncoding_sig <- noncoding_sig[noncoding_sig$Chr == arrayid,]

effect_sizes <- NULL
for(i in 1:nrow(noncoding_sig)){
  ## Chr
  chr <- noncoding_sig$Chr[i]
  ## Gene name
  gene_name <- noncoding_sig$Gene[i]
  ## Coding mask
  category <- noncoding_sig$Category[i]
  
  ### gds file
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  
  a <- Gene_Centric_Noncoding_Burden_Effect_Size_Jake(chr=chr,gene_name=gene_name,category=category ,
                                                 genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
  a <- data.frame(IDK = row.names(a),Gene = gene_name,Chr = chr,Category = category,Number_SNV = a[[4]],cMAC = a[[5]],Burden_Score_Stat = a[[6]],Burden_SE_Score = a[[7]],Burden_pvalue = a[[8]],Burden_Est = a[[9]], Burden_SE_Est = a[[10]])
  effect_sizes <- rbind(effect_sizes,a)
  seqClose(genofile) 
}

noncoding_sig <- noncoding_sig[,-c(1)]
effect_sizes <- effect_sizes[,-c(1)]

noncoding_sig <- inner_join(noncoding_sig,effect_sizes)

write.csv(noncoding_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_noncoding_sigchr",arrayid,".csv"))

