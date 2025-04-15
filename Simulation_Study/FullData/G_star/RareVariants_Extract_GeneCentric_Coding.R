rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricCoding/chr22.Rdata")
coding_sig <- NULL
for(i in 1:length(results_coding)){
  coding_sig <- rbind(coding_sig,unlist(results_coding[[i]][,1:3]))
}

coding_sig <- as.data.frame(coding_sig)
coding_sig$Chr <- as.numeric(coding_sig$Chr)

Gene_Centric_Coding_RareVariants <- function(chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous"),
                                             genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                             QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                             Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){
  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  genes <- genes_info[genes_info[,2]==chr,]  
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  ## get SNV id, position, REF, ALT (whole genome)
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
  
  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  
  rm(filter)
  gc()
  
  ### Gene
  kk <- which(genes[,1]==gene_name)
  
  sub_start_loc <- genes[kk,3]
  sub_end_loc <- genes[kk,4]
  
  is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in]
  
  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
  
  ## plof
  ## Gencode_Exonic
  GENCODE.EXONIC.Category  <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
  
  variant.id.gene <- seqGetData(genofile, "variant.id")  
  
  if(category == "plof"){
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
    variant.id.gene <- variant.id.gene[lof.in.plof]
  }else if(category == "plof_ds"){
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene <- variant.id.gene[lof.in.plof]
  }else if(category == "missense"){
    lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.missense]
  }else if(category == "disruptive_missense"){
    lof.in.ds <- ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene <- variant.id.gene[lof.in.ds]
  }else{
    lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.synonymous]
  }
  
  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
  
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
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  genotype <- Geno
  
  if(dim(genotype)[2] == 1){
    return(matrix(0,nrow = dim(genotype)[1],ncol = 1))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]
  G <- Geno_rare
  rm(Geno_rare)
  gc()
  
  if(is.null(dim(G))){
    G <- matrix(G,ncol = 1)
  }
  
  seqResetFilter(genofile)
  
  return(G)
}

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_Annotation_name_catalog.Rdata"))

obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

rarevariants_gene_centric_coding <- list()

for(i in 1:nrow(coding_sig)){
  chr <- coding_sig$Chr[i]
  ## Gene name
  gene_name <- coding_sig$Gene[i]
  ## Coding mask
  category <- coding_sig$Category[i]
  
  rarevariants_gene_centric_coding[[i]] <- Gene_Centric_Coding_RareVariants(chr=chr,gene_name=gene_name,category=category ,
                                                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE) 
}

seqClose(genofile)

save(rarevariants_gene_centric_coding,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Coding_RareVariants.RData")