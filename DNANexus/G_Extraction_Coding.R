rm(list = ls())
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(data.table)

# for array in {1..22};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/CommonVariant_PRS/G_Extraction_Coding.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/CommonVariant_PRS/G_Extraction_Coding.sh -icmd="bash G_Extraction_Coding.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/ --priority low --instance-type mem3_ssd1_v2_x4
# done

Gene_Centric_Coding_G_Star <- function(chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous"),
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
  
  C <- G%*%matrix(1,nrow=ncol(G),ncol = 1)
  
  seqResetFilter(genofile)
  
  return(C)
}



chr <- as.numeric(commandArgs(TRUE)[1])

Train_PVals_All <- read.csv(paste0("coding_sig.csv"))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$Chr == chr,]

## agds dir

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  ## Null Model
  obj_nullmodel_train <- get(load(paste0(trait,"_Train_Null_Model.RData")))
  obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
  obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))
  
  obj_nullmodel <- obj_nullmodel_train
  obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)
  
  ## Parameter
  QC_label <- "annotation/info/QC_label"
  geno_missing_imputation <- "mean"
  variant_type <- "SNV"	
  
  ## Annotation_dir
  Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
  ## Annotation channel
  Annotation_name_catalog <- get(load("Annotation_name_catalog.Rdata"))
  
  G_star_gene_centric_coding <- list()
  
  gds.path <- list.files()[str_detect(list.files(),".gds")]
  genofile <- seqOpen(gds.path)
  
  for(i in 1:nrow(Train_PVals_All)){
    ## Chr
    chr <- Train_PVals_All$Chr[i]
    ## Gene name
    gene_name <- Train_PVals_All$Gene[i]
    ## Coding mask
    category <- Train_PVals_All$Category[i]
    
    G_star_gene_centric_coding[[i]] <- Gene_Centric_Coding_G_Star(chr=chr,gene_name=gene_name,category=category ,
                                                                  genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE) 
  } 
  
  seqClose(genofile)
  
  G_star_gene_centric_coding <- do.call(cbind,G_star_gene_centric_coding)
  
  fwrite(G_star_gene_centric_coding,file = paste0(trait,"_G_Star_Coding_Chr",chr,".csv"))
}



