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

## source code
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/Burden_Effect_Size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/disruptive_missense_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/missense_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/plof_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/plof_ds_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/synonymous_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/Gene_Centric_Coding_Burden_Effect_Size.R")

###########################################################
#           User Input
###########################################################
## agds dir
agds_dir <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/obj.STAAR.CAD.freeze9.noage.20230620.Rdata"))

## Parameter
QC_label <- "annotation/filter"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/Annotation_name_catalog.Rdata"))
	
###########################################################
#      11, ARFIP2, plof_ds
###########################################################

## Chr
chr <- 11
## Gene name
gene_name <- "ARFIP2"
## Coding mask
category <- "plof_ds"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Coding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
seqClose(genofile)


###########################################################
#      1, SHC1, plof
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "SHC1"
## Coding mask
category <- "plof"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Coding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
seqClose(genofile)


###########################################################
#      5, CDX1, synonymous
###########################################################

## Chr
chr <- 5
## Gene name
gene_name <- "CDX1"
## Coding mask
category <- "synonymous"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Coding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
seqClose(genofile)


###########################################################
#      2, HPCAL1, missense
###########################################################

## Chr
chr <- 2
## Gene name
gene_name <- "HPCAL1"
## Coding mask
category <- "missense"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Coding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
seqClose(genofile)


###########################################################
#      15, CAPN3, disruptive_missense
###########################################################

## Chr
chr <- 15
## Gene name
gene_name <- "CAPN3"
## Coding mask
category <- "disruptive_missense"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Coding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
seqClose(genofile)






		
