#####################################################################
# Gene-centric analysis for noncoding rare variants of protein-coding 
# genes using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
#####################################################################
rm(list=ls())
gc()

# for trait in 1 2 3 4 5 6;
# do
# for array in {1..379};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/RareVariant_Analysis_Noncoding.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/RareVariant_Analysis_Noncoding.sh -icmd="bash RareVariant_Analysis_Noncoding.sh ${trait} ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2
# done
# done

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "BMI"
}else if(trait == 2){
  trait <- "TC"
}else if(trait == 3){
  trait <- "HDL"
}else if(trait == 4){
  trait <- "LDL"
}else if(trait == 5){
  trait <- "logTG"
}else{
  trait <- "Height"
}

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(readr)
library(dplyr)
library(stringr)

###########################################################
#           User Input
###########################################################
## aGDS directory
## Null model
obj_nullmodel <- get(load(paste0(trait,"_Train_Null_Model.RData")))

## Parameter
QC_label <- "annotation/info/QC_label2"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output file name
output_file_name <- paste0(trait,"_UKBB_WGS_Noncoding_Train")
## input array id from batch file (Harvard FAS RC cluster)
arrayid <- as.numeric(commandArgs(TRUE)[2])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50 
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)

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

## exclude large noncoding masks
jobid_exclude <- c(21,39,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
sub_seq_id_exclude <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)

for(i in 1:length(jobid_exclude)){
  if(arrayid==jobid_exclude[i]){
    sub_seq_id <- setdiff(sub_seq_id,sub_seq_id_exclude[i])
  }
}

## aGDS file
agds.path <- list.files()[str_detect(list.files(),".gds")]
genofile <- seqOpen(agds.path)

genes <- genes_info

results_noncoding <- c()
for(kk in sub_seq_id){
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_noncoding <- append(results_noncoding,results)
}

save(results_noncoding,file=paste0(output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

system(paste0("rm Annotation_name_catalog.csv"))
system(paste0("rm ",agds.path))
system(paste0("rm ",paste0(trait,"_Train_Null_Model.RData")))
