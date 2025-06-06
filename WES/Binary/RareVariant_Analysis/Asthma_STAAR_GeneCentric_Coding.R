#####################################################################
# Gene-centric analysis for coding rare variants using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
#####################################################################

rm(list=ls())
gc()

trait <- "Asthma"

## load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################

## job nums
jobs_num <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Train_Null_Model.RData")))

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/Annotation_name_catalog.Rdata"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/"
## output file name
output_file_name <- paste0(trait,"_UKBB_WES_Coding_Train")
## input array id from batch file (Harvard FAS cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])
# arrayid <- 1

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50 
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
group.num <- group.num.allchr[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num){ 
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}	

### exclude large genes
if(arrayid==57){
  sub_seq_id <- setdiff(sub_seq_id,840)
}

if(arrayid==112){
  sub_seq_id <- setdiff(sub_seq_id,c(543,544))
}

if(arrayid==113){
  sub_seq_id <- setdiff(sub_seq_id,c(575,576,577,578,579,580))
}

### gds file
gds.dir <- "/data/williamsjacr/UKB_WES_Full_Processed_Data/gds/"
gds.path <- paste0(gds.dir,"chr",chr,".gds")
genofile <- seqOpen(gds.path)

genes <- genes_info

results_coding <- c()

for(kk in sub_seq_id){
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_coding <- append(results_coding,results)
}

save(results_coding,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)
