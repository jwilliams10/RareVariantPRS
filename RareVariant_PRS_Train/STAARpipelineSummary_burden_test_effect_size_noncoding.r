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

## source code
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/Burden_Effect_Size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/UTR_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/upstream_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/downstream_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/promoter_CAGE_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/promoter_DHS_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/enhancer_CAGE_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/enhancer_DHS_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/ncRNA_burden_effect_size.R")
source("/n/holystore01/LABS/xlin/Lab/xihao_zilin/TOPMed_CAD_F9_noage_noPROMIS/burden_effect_size_test/Gene_Centric_Noncoding_Burden_Effect_Size.R")

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
#      18, GREB1L, downstream
###########################################################

## Chr
chr <- 18
## Gene name
gene_name <- "GREB1L"
## Coding mask
category <- "downstream"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)



Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)

## Burden p: 0.06013495

seqClose(genofile)


###########################################################
#      3, GRIP2, enhancer_DHS
###########################################################

## Chr
chr <- 3
## Gene name
gene_name <- "GRIP2"
## Coding mask
category <- "enhancer_DHS"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
## Burden p: 0.3828957

seqClose(genofile)


###########################################################
#      1, OR2M4, UTR
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "OR2M4"
## Coding mask
category <- "UTR"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)
								
## Burden p: 0.6120259

seqClose(genofile)


###########################################################
#      1, OR2M4, upstream
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "OR2M4"
## Coding mask
category <- "upstream"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)

## Burden p: 0.9195464									

seqClose(genofile)

	
###########################################################
#      1, RABGGTB, enhancer_CAGE
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "RABGGTB"
## Coding mask
category <- "enhancer_CAGE"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)

## Burden p: 0.7001279									
								
seqClose(genofile)

	
###########################################################
#      1, RABGGTB, promoter_CAGE
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "RABGGTB"
## Coding mask
category <- "promoter_CAGE"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)

## Burden p: 0.9214745									
								
seqClose(genofile)		
		

###########################################################
#      1, RABGGTB, promoter_DHS
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "RABGGTB"
## Coding mask
category <- "promoter_DHS"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)

## Burden p: 0.9198949									
								
seqClose(genofile)		
		

###########################################################
#      1, MIR6859-1, ncRNA
###########################################################

## Chr
chr <- 1
## Gene name
gene_name <- "MIR6859-1"
## Coding mask
category <- "ncRNA"

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

results <- Gene_Centric_Noncoding_Burden_Effect_Size(chr=chr,gene_name=gene_name,category=category ,
                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)

results

## Burden p: 0.5877116		

seqClose(genofile)



