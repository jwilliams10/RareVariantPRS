
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
library(readr)
library(dplyr)

trait <- "BMI"

## source code
Burden_Effect_Size <- function(genotype,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2){
  
  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]
  
  rm(genotype)
  gc()
  
  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()
    
    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        
        ## dense GRM
        P <- obj_nullmodel$P
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        
        Cov <- t(G)%*%P%*%G
      }else{
        ## Sparse GRM 
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov
        
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        
        tSigma_iX_G = t(Sigma_iX)%*%G
        Cov = t(Sigma_i%*%G)%*%G - t(tSigma_iX_G)%*%cov%*%tSigma_iX_G
        
      }
    }else{
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        tX_G = t(X)%*%diag(working)%*%G
        Cov = t(G)%*%diag(working)%*%G - t(tX_G)%*%solve(t(X)%*%diag(working)%*%X)%*%tX_G
        
      }else if(obj_nullmodel$family[1] == "gaussian"){
        tX_G = t(X)%*%G
        Cov = t(G)%*%G - t(tX_G)%*%solve(t(X)%*%X)%*%tX_G
      }
      
      Cov = Cov*sigma^2
      
      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      
    }
    
    Score <- t(G)%*%residuals.phenotype 
    
    Burden_Score <- sum(Score)
    Burden_Var <- sum(Cov)
    
    Burden_Effect_Size <- Burden_Score/Burden_Var
    Burden_test_chisq <- Burden_Score^2/Burden_Var
    Burden_pvalue <- pchisq(Burden_test_chisq,1,lower=FALSE)
    
    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G)
    
    return(data.frame(num_variant = num_variant,
                      cMAC = cMAC,
                      RV_label = RV_label,
                      Burden_Score_Stat = Burden_Score, 
                      Burden_SE_Score = sqrt(Burden_Var),
                      Burden_pvalue = Burden_pvalue,
                      Burden_Est = Burden_Effect_Size,
                      Burden_SE_Est = 1/sqrt(Burden_Var)))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }
  
}




Gene_Centric_Coding_Burden_Effect_Size_Jake <- function(category_filter,sub_start_loc,sub_end_loc,variant.id,position,SNVlist,id.genotype.match,
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
  
  burden_eff <- NULL
  try(burden_eff <- Burden_Effect_Size(genotype=Geno,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  
  seqResetFilter(genofile)
  return(burden_eff)
}




###########################################################
#           User Input
###########################################################

### Significant Results 
coding_sig <- read_csv(paste0(trait,"_coding_sig.csv"))
colnames(coding_sig) <- c("IDK","Gene","Chr","Category","Number_SNV","SKAT_1_25","Burden_1_1","ACAT_V_1_25","STAAR_O")

## agds dir
agds_dir <- get(load("agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load(paste0(trait,"_Train_Null_Model.RData")))

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("Annotation_name_catalog.Rdata"))

chunks <- list()

set <- c(0,22,35,47,55,64,75,85,92,100,108,122,133,136,143,149,158,170,173,189,194,196,200)

for(i in 1:22){
  indexes <- which(coding_sig$Chr == i)
  a <- cut(seq_along(indexes), set[i + 1] - set[i], labels = FALSE)
  chunks <- c(chunks,split(indexes,a))
}

arrayid <- 20

coding_sig <- coding_sig[chunks[[arrayid]],]


##############################
## Preprocessing
##############################

chr <- unique(coding_sig$Chr)

gds.path <- agds_dir[chr]
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
for(i in 1:nrow(coding_sig)){
  kk_vector[i] <- which(genes[,1]==coding_sig$Gene[i])
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

unique_ids <- paste0(coding_sig$Gene,"__",coding_sig$Category,"__",sub_start_loc_vec,"__",sub_end_loc_vec)
x <- unique_ids[1]

coding_effectsizes_parlapply <- function(x,
                                         category_filter_plof,category_filter_plof_ds,category_filter_missense,category_filter_disruptive_missense,category_filter_synonymous,
                                         variant.id,id.genotype.match,genofile,obj_nullmodel,rare_maf_cutoff,rv_num_cutoff,geno_missing_imputation){
  tmp <- unlist(strsplit(x,"__"))
  gene <- tmp[1]
  category <- tmp[2]
  sub_start_loc <- as.numeric(tmp[3])
  sub_end_loc <- as.numeric(tmp[4])
  
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
  
  a <- Gene_Centric_Coding_Burden_Effect_Size_Jake(sub_start_loc = sub_start_loc,sub_end_loc = sub_end_loc,variant.id = variant.id,id.genotype.match = id.genotype.match,
                                                   SNVlist = SNVlist,position = position, category_filter = category_filter, genofile = genofile,obj_nullmodel = obj_nullmodel,
                                                   rare_maf_cutoff=0.01,rv_num_cutoff=2,geno_missing_imputation=geno_missing_imputation,silent=FALSE)
  
  if(!is.null(a)){
    a <- data.frame(Gene = gene, Category = category,Number_SNV = a$num_variant,cMAC = a$cMAC,Burden_Score_Stat = a$Burden_Score_Stat,
                    Burden_SE_Score = a$Burden_SE_Score,Burden_pvalue = a$Burden_pvalue,Burden_Est = a$Burden_Est, Burden_SE_Est = a$Burden_SE_Est)
    a <- unique(a)
  }
  return(a)
} 

a <- lapply(unique_ids,coding_effectsizes_parlapply,
            category_filter_plof = category_filter_plof,category_filter_plof_ds = category_filter_plof_ds,
            category_filter_missense = category_filter_missense,category_filter_disruptive_missense = category_filter_disruptive_missense,
            category_filter_synonymous = category_filter_synonymous, variant.id = variant.id,id.genotype.match = id.genotype.match,genofile = genofile,
            obj_nullmodel = obj_nullmodel,rare_maf_cutoff = rare_maf_cutoff,rv_num_cutoff = rv_num_cutoff,geno_missing_imputation = geno_missing_imputation)

seqClose(genofile)

effect_sizes <- data.table::rbindlist(a)
effect_sizes$Chr <- unique(coding_sig$Chr)

coding_sig <- coding_sig[,-c(1)]

coding_sig <- inner_join(coding_sig,effect_sizes)
