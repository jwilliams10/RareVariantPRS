rm(list=ls())
gc()

# for trait in 1 2 3 4 5 6;
# do
# for array in {1..379};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/RareVariant_Analysis_Noncoding_Preload.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/RareVariant_Analysis_Noncoding_Preload.sh -icmd="bash RareVariant_Analysis_Noncoding_Preload.sh ${trait} ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2
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

if(chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
noncoding_longmasks <- c("ARID4B","RERE","CRIM1","EML4","EML6","GLI2","ITGB6","FNDC3B","AFAP1","TRIO",
                         "NR3C1","SGK1","PRDM1","TNS3","GNA12","HOXA2","PLEC","ASAP1","CSGALNACT1","RAPGEF1",
                         "KLF6","DNAJB12","ABLIM1","ETV6","IGF1R","RMI2","CMIP","HOXB2","SEPT9","SPECC1",
                         "CIRBP","ADNP","CDC42EP3","EPAS1","MGAT5")
genes_info_chr <- genes_info_chr[!genes_info_chr[,1] %in% noncoding_longmasks,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num){
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}

Gene_Centric_Noncoding_preload <- function(chr,gene_name,category=c("all_categories","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                           genofile,obj_nullmodel,
                                           dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                           dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                           dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                           dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                           QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                           Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                           SPA_p_filter=TRUE,p_filter_cutoff=0.05,silent=FALSE){
  
  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  if(category=="all_categories"){
    results <- noncoding_preload(chr,gene_name,genofile,obj_nullmodel,
                                 dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                 dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                 dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                 dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="downstream"){
    results <- downstream(chr,gene_name,genofile,obj_nullmodel,
                          rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                          SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="upstream"){
    results <- upstream(chr,gene_name,genofile,obj_nullmodel,
                        rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                        SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="UTR"){
    results <- UTR(chr,gene_name,genofile,obj_nullmodel,
                   rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                   QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                   SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="promoter_CAGE"){
    results <- promoter_CAGE_preload(chr,gene_name,genofile,obj_nullmodel,
                                     dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                     rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                     SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="promoter_DHS"){
    results <- promoter_DHS_preload(chr,gene_name,genofile,obj_nullmodel,
                                    dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                    rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="enhancer_CAGE"){
    results <- enhancer_CAGE_preload(chr,gene_name,genofile,obj_nullmodel,
                                     dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                     rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                     SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  if(category=="enhancer_DHS"){
    results <- enhancer_DHS_preload(chr,gene_name,genofile,obj_nullmodel,
                                    dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                    rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }
  
  return(results)
}

noncoding_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                              dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                              dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                              dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                              dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                              rare_maf_cutoff=0.01,rv_num_cutoff=2,
                              QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                              Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                              Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                              SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){
  
  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno
  
  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA)){
    use_SPA <- obj_nullmodel$use_SPA
  }else{
    use_SPA <- FALSE
  }
  
  #####################################
  #   Gene Info
  ## get SNV id
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant"){
    SNVlist <- filter == "PASS"
  }
  
  if(variant_type=="SNV"){
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  
  if(variant_type=="Indel"){
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  
  variant.id <- seqGetData(genofile, "variant.id")
  
  rm(filter)
  gc()
  
  ########################################
  #   Downstream
  ## downstream SNVs
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- (GENCODE.Category=="downstream")&(SNVlist)
  variant.id.downstream <- variant.id[is.in]
  
  seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)
  
  rm(variant.id.downstream)
  gc()
  
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
  
  rm(GENCODE.Info)
  gc()
  
  rm(variant_gene_num)
  gc()
  
  Gene <- as.character(unlist(GENCODE.Info.split))
  
  rm(GENCODE.Info.split)
  gc()
  
  seqResetFilter(genofile)
  
  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_downstream <- c()
  if(inherits(pvalues, "list")){
    results_temp <- rep(NA,4)
    results_temp[3] <- "downstream"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_downstream <- rbind(results_downstream,results_temp)
  }
  
  if(!is.null(results_downstream)){
    if(!use_SPA){
      colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
      colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
      colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_downstream)[dim(results_downstream)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  ########################################
  #   Upstream
  
  is.in <- (GENCODE.Category=="upstream")&(SNVlist)
  variant.id.upstream <- variant.id[is.in]
  
  seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)
  
  rm(variant.id.upstream)
  gc()
  
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
  
  rm(GENCODE.Info)
  gc()
  
  rm(variant_gene_num)
  gc()
  
  Gene <- as.character(unlist(GENCODE.Info.split))
  
  rm(GENCODE.Info.split)
  gc()
  
  seqResetFilter(genofile)
  
  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_upstream <- c()
  if(inherits(pvalues, "list")){
    results_temp <- rep(NA,4)
    results_temp[3] <- "upstream"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_upstream <- rbind(results_upstream,results_temp)
  }
  
  if(!is.null(results_upstream)){
    if(!use_SPA){
      colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
      colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
      colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_upstream)[dim(results_upstream)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  ########################################################
  #                UTR
  
  is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
  variant.id.UTR <- variant.id[is.in]
  
  rm(GENCODE.Category)
  gc()
  
  seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)
  
  rm(variant.id.UTR)
  gc()
  
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
  
  rm(GENCODE.Info)
  gc()
  
  # Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[seq(1,length(z),2)]))
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))
  
  rm(GENCODE.Info.split)
  gc()
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  
  seqResetFilter(genofile)
  
  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_UTR <- c()
  if(inherits(pvalues, "list")){
    results_temp <- rep(NA,4)
    results_temp[3] <- "UTR"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_UTR <- rbind(results_UTR,results_temp)
  }
  
  if(!is.null(results_UTR)){
    if(!use_SPA){
      colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
      colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
      colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_UTR)[dim(results_UTR)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  #############################################
  #   Promoter-CAGE
  
  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromCAGE[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_promoter_CAGE <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfPromCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
  }
  
  if(!is.null(results_promoter_CAGE)){
    if(!use_SPA){
      colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_CAGE)[dim(results_promoter_CAGE)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  ##################################################
  #       Promoter-DHS
  
  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromrOCRs[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_promoter_DHS <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfPromrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_promoter_DHS <- rbind(results_promoter_DHS ,results_temp)
  }
  
  if(!is.null(results_promoter_DHS)){
    if(!use_SPA){
      colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_promoter_DHS)[dim(results_promoter_DHS)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  ###########################################
  #        Enhancer-CAGE
  
  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerCAGE[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_enhancer_CAGE <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfHancerCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
  }
  
  if(!is.null(results_enhancer_CAGE)){
    if(!use_SPA){
      colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_CAGE)[dim(results_enhancer_CAGE)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  ##################################################
  #       Enhancer-DHS
  
  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerrOCRs[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results_enhancer_DHS <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfHancerrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
  }
  
  if(!is.null(results_enhancer_DHS)){
    if(!use_SPA){
      colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
      colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results_enhancer_DHS)[dim(results_enhancer_DHS)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  ############################################
  #           results
  
  results_noncoding <- list(upstream=results_upstream,downstream=results_downstream,UTR=results_UTR,
                            promoter_CAGE=results_promoter_CAGE,promoter_DHS=results_promoter_DHS,
                            enhancer_CAGE=results_enhancer_CAGE,enhancer_DHS=results_enhancer_DHS)
  
  return(results_noncoding)
  
}

promoter_CAGE_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                  dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                  SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){
  
  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno
  
  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA)){
    use_SPA <- obj_nullmodel$use_SPA
  }else{
    use_SPA <- FALSE
  }
  
  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromCAGE[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfPromCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    
    results <- rbind(results,results_temp)
  }
  
  if(!is.null(results)){
    if(!use_SPA){
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  return(results)
}

promoter_DHS_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                 dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                 Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                 SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){
  
  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno
  
  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA)){
    use_SPA <- obj_nullmodel$use_SPA
  }else{
    use_SPA <- FALSE
  }
  
  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromrOCRs[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfPromrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "promoter_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results <- rbind(results,results_temp)
  }
  
  
  if(!is.null(results)){
    if(!use_SPA){
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  return(results)
}

enhancer_CAGE_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                  dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                  SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){
  
  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno
  
  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA)){
    use_SPA <- obj_nullmodel$use_SPA
  }else{
    use_SPA <- FALSE
  }
  
  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerCAGE[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfHancerCAGEVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_CAGE"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results <- rbind(results,results_temp)
  }
  
  if(!is.null(results)){
    if(!use_SPA){
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  return(results)
}

enhancer_DHS_preload <- function(chr,gene_name,genofile,obj_nullmodel,
                                 dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                 Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                 SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){
  
  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno
  
  ## SPA status
  if(!is.null(obj_nullmodel$use_SPA)){
    use_SPA <- obj_nullmodel$use_SPA
  }else{
    use_SPA <- FALSE
  }
  
  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerrOCRs[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV"){
    if(Use_annotation_weights){
      for(k in 1:length(Annotation_name)){
        if(Annotation_name[k]%in%Annotation_name_catalog$name){
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD"){
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity"){
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  pvalues <- 0
  if(n_pheno == 1){
    if(!use_SPA){
      try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
    }else{
      try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
    }
  }else{
    try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  }
  
  results <- c()
  if(inherits(pvalues, "list")){
    results_temp <- dfHancerrOCRsVarGene.SNV[1,1:4]
    results_temp[3] <- "enhancer_DHS"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    
    if(!use_SPA){
      results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                        pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
    }else{
      results_temp <- c(results_temp,pvalues$cMAC,
                        pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
    }
    
    results <- rbind(results,results_temp)
  }
  
  if(!is.null(results)){
    if(!use_SPA){
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
    }else{
      colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
      colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }
  
  seqResetFilter(genofile)
  
  return(results)
}

gene_centric_noncoding_dnanexus <- function(genes_info_chr,gene_name,chr,genofile,obj_nullmodel,
                                            dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                            dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                            dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                            dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                            rare_maf_cutoff,QC_label,variant_type,geno_missing_imputation,
                                            Annotation_dir,Annotation_name_catalog,
                                            Use_annotation_weights,Annotation_name,silent){
  results <- try(Gene_Centric_Noncoding_preload(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                                dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                                dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                                dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                                rare_maf_cutoff=rare_maf_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent),silent=TRUE)
  return(results)
}

#########################################################
#             Promoter_CAGE
#########################################################

agds.path <- list.files()[str_detect(list.files(),".gds")]

genofile <- seqOpen(agds.path)

varid <- seqGetData(genofile, "variant.id")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

#Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
CAGEBvt <- CAGEAnno!=""
CAGEidx <- which(CAGEBvt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[CAGEidx])
seqSetFilter(genofile,promGobj,intersect=TRUE)
CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
##obtain variants info
CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
CAGEvref <- as.character(seqGetData(genofile,"$ref"))
CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant"){
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV"){
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel"){
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.PromCAGE <- variant.id[SNVlist]

dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

seqResetFilter(genofile)

rm(dfPromCAGEVarGene)
gc()

#########################################################
#             Promoter_DHS
#########################################################

varid <- seqGetData(genofile, "variant.id")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

# Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRsBvt <- rOCRsAnno!=""
rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[rOCRsidx])

seqSetFilter(genofile,promGobj,intersect=TRUE)
rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
## obtain variants info
rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant"){
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV"){
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel"){
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.PromrOCRs <- variant.id[SNVlist]

dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

seqResetFilter(genofile)

rm(dfPromrOCRsVarGene)
gc()

#########################################################
#             Enhancer_CAGE
#########################################################

varid <- seqGetData(genofile, "variant.id")

#Now extract the GeneHancer with CAGE Signal Overlay
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
genehancer <- genehancerAnno!=""

CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
CAGE <- CAGEAnno!=""
CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

# variants that covered by whole GeneHancer without CAGE overlap.
genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
enhancervpos <- as.numeric(seqGetData(genofile,"position"))
enhancervref <- as.character(seqGetData(genofile,"$ref"))
enhancervalt <- as.character(seqGetData(genofile,"$alt"))
dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant"){
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV"){
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel"){
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.HancerCAGE <- variant.id[SNVlist]

dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)

seqResetFilter(genofile)

rm(dfHancerCAGEVarGene)
gc()

#########################################################
#             Enhancer_DHS
#########################################################

varid <- seqGetData(genofile, "variant.id")

#Now extract the GeneHancer with rOCRs Signal Overlay
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
genehancer <- genehancerAnno!=""

rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRs <- rOCRsAnno!=""
rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
# variants that covered by whole GeneHancer without rOCRs overlap.

genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
enhancervpos <- as.numeric(seqGetData(genofile,"position"))
enhancervref <- as.character(seqGetData(genofile,"$ref"))
enhancervalt <- as.character(seqGetData(genofile,"$alt"))
dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

rm(varid)
gc()

## get SNV id
filter <- seqGetData(genofile, QC_label)
if(variant_type=="variant"){
  SNVlist <- filter == "PASS"
}

if(variant_type=="SNV"){
  SNVlist <- (filter == "PASS") & isSNV(genofile)
}

if(variant_type=="Indel"){
  SNVlist <- (filter == "PASS") & (!isSNV(genofile))
}

variant.id <- seqGetData(genofile, "variant.id")
variant.id.SNV.HancerrOCRs <- variant.id[SNVlist]

dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)

seqResetFilter(genofile)

rm(dfHancerrOCRsVarGene)
gc()

results_noncoding <- c()
for(kk in sub_seq_id){
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- gene_centric_noncoding_dnanexus(genes_info_chr=genes_info_chr,gene_name = gene_name,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                             dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                             dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                             dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                             dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                             rare_maf_cutoff=0.01,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=TRUE)
  
  results_noncoding <- append(results_noncoding,results)
}

seqClose(genofile)

rm(list=setdiff(ls(), c("results_noncoding", "output_file_name"))); gc()
save(results_noncoding, file = paste0(output_file_name,"_",arrayid,".Rdata"))

system(paste0("rm Annotation_name_catalog.csv"))
system(paste0("rm ",agds.path))
system(paste0("rm ",paste0(trait,"_Train_Null_Model.RData")))