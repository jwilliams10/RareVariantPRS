rm(list=ls())
list.files()
gc()

# for trait in 1 2 3 4 5 6;
# do
# for array in {1..1808};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/GeneCentric_Coding_PRS.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/GeneCentric_Coding_PRS.sh -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/coding_chr.csv -icmd="bash GeneCentric_Coding_PRS.sh ${trait} ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_Coding_PRS/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2
# done
# done

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
library(stringr)

Burden_PRS <- function(genotype,rare_maf_cutoff=0.01,rv_num_cutoff=2,BETA = NULL){
  
  # if(class(genotype) != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
  #   stop("genotype is not a matrix!")
  # }
  
  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  
  gc()
  
  if(sum(RV_label) == 0){
    PRS <- matrix(0,nrow = nrow(genotype$Geno),ncol = 1)
    return(PRS)
  }else if(sum(RV_label) == 1){
    Geno_rare <- genotype$Geno[,RV_label,drop = FALSE]
    G <- Geno_rare
    rm(Geno_rare)
    gc()
    
    PRS <- G*BETA
    return(PRS)
  }else{
    Geno_rare <- genotype$Geno[,RV_label]
    G <- Geno_rare
    rm(Geno_rare)
    gc()
    
    C <- G%*%matrix(1,nrow=ncol(G),ncol = 1)
    PRS <- C*BETA
    return(PRS)
  }
  
}

Gene_Centric_Coding_Burden_PRS <- function(category_filter,BETA,sub_start_loc,sub_end_loc,variant.id,position,SNVlist,id.genotype.match,
                                                        genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,geno_missing_imputation=c("mean","minor"),silent=FALSE){
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in & category_filter]
  
  seqSetFilter(genofile,variant.id=variant.id.gene)
  
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
  
  Burden_PRS <- Burden_PRS(genotype=Geno,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,BETA = BETA)
  
  Burden_PRS <- data.frame(ID = phenotype.id.merge$phenotype.id,PRS = Burden_PRS)
  
  seqResetFilter(genofile)
  return(Burden_PRS)
}

###########################################################
#           User Input
###########################################################

trait <- gsub("_coding_sig_es.csv","",list.files()[str_detect(list.files(),"_coding_sig_es.csv")])

### Significant Results 
Train_Effect_Sizes_All <- read_csv(paste0(trait,"_coding_sig_es.csv"))

## Null Model
obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_tune
obj_nullmodel$id_include <- c(obj_nullmodel$id_include,obj_nullmodel_validation$id_include)

## Parameter
QC_label <- "annotation/info/QC_label2"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

PRS <- NULL

arrayid <- as.numeric(commandArgs(TRUE)[1])
arrayid_original <- arrayid

if(arrayid>904){
  Burden <- 1
  arrayid <- arrayid - 904
}else{
  Burden <- 0
}

sets <- c(0,22,44,66,88,110,132,154,176,198,220,242,264,304,504,904)
splits <- diff(sets)

if(arrayid <= 22){
  threshold <- 1
  jobid <- arrayid - 0
}else if(arrayid <= 44){
  threshold <- 2
  jobid <- arrayid - 22
}else if(arrayid <= 66){
  threshold <- 3
  jobid <- arrayid - 44
}else if(arrayid <= 88){
  threshold <- 4
  jobid <- arrayid - 66
}else if(arrayid <= 110){
  threshold <- 5
  jobid <- arrayid - 88
}else if(arrayid <= 132){
  threshold <- 6
  jobid <- arrayid - 110
}else if(arrayid <= 154){
  threshold <- 7
  jobid <- arrayid - 132
}else if(arrayid <= 176){
  threshold <- 8
  jobid <- arrayid - 154
}else if(arrayid <= 198){
  threshold <- 9
  jobid <- arrayid - 176
}else if(arrayid <= 220){
  threshold <- 10
  jobid <- arrayid - 198
}else if(arrayid <= 242){
  threshold <- 11
  jobid <- arrayid - 220
}else if(arrayid <= 264){
  threshold <- 12
  jobid <- arrayid - 242
}else if(arrayid <= 304){
  threshold <- 13
  jobid <- arrayid - 264
  set <- c(0,4,7,9,11,13,15,17,18,20,22,24,27,28,29,30,32,34,35,37,38,39,40)
}else if(arrayid <= 504){
  threshold <- 14
  jobid <- arrayid - 304
  set <- c(0,22,35,47,55,64,75,85,92,100,108,122,133,136,143,149,158,170,173,189,194,196,200)
}else{
  threshold <- 15
  jobid <- arrayid - 504
  set <- c(0,44,70,94,110,128,150,170,184,200,216,244,266,272,286,298,316,340,346,378,388,392,400)
}

if(Burden == 0){
  Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$STAAR_O <= thresholds[threshold],]
}else{
  Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$Burden_1_1 <= thresholds[threshold],]
}

Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$Burden_Est != 0,]

if(nrow(Train_Effect_Sizes_All) == 0){
  PRS <- data.frame(ID = 1:length(obj_nullmodel$id_include),PRS = 0)
}else{
  chunks <- NULL
  if(splits[threshold] != 22){
    for(i in 1:22){
      indexes <- which(Train_Effect_Sizes_All$Chr == i)
      if(set[i + 1] - set[i] != 1){
        a <- cut(seq_along(indexes), set[i + 1] - set[i], labels = FALSE)
        chunks <- c(chunks,split(indexes,a))
      }else{
        chunks <- c(chunks,list(indexes))
      }
    }
    Train_Effect_Sizes_All <- Train_Effect_Sizes_All[chunks[[jobid]],]
  }else{
    Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$Chr == jobid,]
  }
  
  if(nrow(Train_Effect_Sizes_All) == 0){
    PRS <- data.frame(ID = 1:length(obj_nullmodel$id_include),PRS = 0)
  }else{
    ##############################
    ## Preprocessing
    ##############################
    
    chr <- unique(Train_Effect_Sizes_All$Chr)
    
    print(chr)
    
    gds.path <- list.files()[str_detect(list.files(),".gds")]
    genofile <- seqOpen(gds.path)
    
    phenotype.id <- as.character(obj_nullmodel$id_include)
    
    genes <- genes_info[genes_info[,2]==chr,]  
    
    ## get SNV id, position, REF, ALT (whole genome)
    filter <- seqGetData(genofile, QC_label)
    SNVlist <- (filter == "PASS") & isSNV(genofile)
    
    position <- as.numeric(seqGetData(genofile, "position"))
    variant.id <- seqGetData(genofile, "variant.id")
    
    rm(filter)
    gc()
    
    kk_vector <- vector()
    for(i in 1:nrow(Train_Effect_Sizes_All)){
      kk_vector[i] <- which(genes[,1]==Train_Effect_Sizes_All$Gene[i])
    }
    
    sub_start_loc_vec <- genes[kk_vector,3]
    sub_end_loc_vec <- genes[kk_vector,4]
    
    ## Gencode_Exonic
    GENCODE.EXONIC.Category  <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
    ## Gencode
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    ## Meta.SVM.Pred
    MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
    
    category_filter_plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
    category_filter_plof_ds <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    category_filter_missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
    category_filter_disruptive_missense <- ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    category_filter_synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
    
    id.genotype <- seqGetData(genofile,"sample.id")
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    unique_ids <- paste0(Train_Effect_Sizes_All$Gene,"__",Train_Effect_Sizes_All$Category,"__",sub_start_loc_vec,"__",sub_end_loc_vec,"__",Train_Effect_Sizes_All$Burden_Est)
    
    GeneCentricCoding_effectsizes_parlapply <- function(x,
                                                        category_filter_plof,category_filter_plof_ds,category_filter_missense,category_filter_disruptive_missense,category_filter_synonymous,
                                                        variant.id,id.genotype.match,genofile,obj_nullmodel,geno_missing_imputation){
      
      tmp <- unlist(strsplit(x,"__"))
      gene <- tmp[1]
      category <- tmp[2]
      sub_start_loc <- as.numeric(tmp[3])
      sub_end_loc <- as.numeric(tmp[4])
      BETA <- as.numeric(tmp[5])
      
      if(category == "plof"){
        category_filter <- category_filter_plof
      }else if(category == "plof_ds"){
        category_filter <- category_filter_plof_ds
      }else if(category == "missense"){
        category_filter <- category_filter_missense
      }else if(category == "disruptive_missense"){
        category_filter <- category_filter_disruptive_missense
      }else{
        category_filter <- category_filter_synonymous
      }
      
      a <- Gene_Centric_Coding_Burden_PRS(sub_start_loc = sub_start_loc,sub_end_loc = sub_end_loc,variant.id = variant.id,id.genotype.match = id.genotype.match,BETA = BETA,
                                                       SNVlist = SNVlist,position = position, category_filter = category_filter, genofile = genofile,obj_nullmodel = obj_nullmodel,
                                                       rare_maf_cutoff=0.01,rv_num_cutoff=2,geno_missing_imputation=geno_missing_imputation,silent=FALSE)
      return(a)
    } 
    
    PRSs <- lapply(unique_ids,GeneCentricCoding_effectsizes_parlapply,
                   category_filter_plof = category_filter_plof,category_filter_plof_ds = category_filter_plof_ds,
                   category_filter_missense = category_filter_missense,category_filter_disruptive_missense = category_filter_disruptive_missense,
                   category_filter_synonymous = category_filter_synonymous, variant.id = variant.id,id.genotype.match = id.genotype.match,genofile = genofile,
                   obj_nullmodel = obj_nullmodel,geno_missing_imputation = geno_missing_imputation)
    
    seqClose(genofile)
    
    PRS_Final <- PRSs[[1]]
    if(length(PRSs) > 1){
      for(i in 2:length(PRSs)){
        PRS_Final[,2] <- PRS_Final[,2] + PRSs[[i]][,2]
      } 
    }
    PRS <- PRS_Final 
  }
}

write.csv(PRS,row.names = FALSE,file = paste0(trait,"_Coding_PRS",arrayid_original,".csv"))

a <- list.files()[list.files() != paste0(trait,"_Coding_PRS",arrayid_original,".csv")]

for(i in 1:length(a)){
  system(paste0("rm ",a[i]))
}
