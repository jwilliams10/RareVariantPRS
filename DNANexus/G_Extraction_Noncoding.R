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
library(stringr)

# for array in {1..22};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/G_Extraction_Noncoding.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/G_Extraction_Noncoding.sh -icmd="bash G_Extraction_Noncoding.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/ --priority low --instance-type mem3_ssd1_v2_x4
# done

Gene_Centric_Noncoding_G_Star <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
                                          genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                          QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                          Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){
  
  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  genes <- genes_info[genes_info[,2]==chr,]
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  
  if(category == "downstream"){
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
    
    ## downstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="downstream")&(SNVlist)
    variant.id.downstream <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)
    rm(variant.id.downstream)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
    
    rm(GENCODE.Info)
    rm(variant_gene_num)
    
    Gene <- as.character(unlist(GENCODE.Info.split))
    
    rm(GENCODE.Info.split)
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category == "upstream"){
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
    
    ## upstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="upstream")&(SNVlist)
    variant.id.upstream <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)
    rm(variant.id.upstream)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
    
    rm(GENCODE.Info)
    rm(variant_gene_num)
    
    Gene <- as.character(unlist(GENCODE.Info.split))
    rm(GENCODE.Info.split)
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category == "UTR"){
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
    
    ## downstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
    variant.id.UTR <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)
    rm(variant.id.UTR)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
    
    rm(GENCODE.Info)
    
    # Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[seq(1,length(z),2)]))
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))
    rm(GENCODE.Info.split)
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category == "promoter_CAGE"){
    ## Promoter
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
    
    print(SNVlist)
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
    dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
    dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
    dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)
    
    seqResetFilter(genofile)
    rm(dfPromCAGEVarGene)
    
    ### Gene
    is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category=="promoter_DHS"){
    ## Promoter
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
    variant.id.SNV <- variant.id[SNVlist]
    
    dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
    dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
    dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
    dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)
    
    seqResetFilter(genofile)
    rm(dfPromrOCRsVarGene)
    
    ### Gene
    is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category=="enhancer_CAGE"){
    ## Enhancer
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
    variant.id.SNV <- variant.id[SNVlist]
    
    dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
    dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
    dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
    dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)
    
    seqResetFilter(genofile)
    
    rm(dfHancerCAGEVarGene)
    
    ### Gene
    is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category=="enhancer_DHS"){
    ## Enhancer
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
    variant.id.SNV <- variant.id[SNVlist]
    
    dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
    dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
    dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
    dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)
    
    seqResetFilter(genofile)
    rm(dfHancerrOCRsVarGene)
    
    ### Gene
    is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else{
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
    
    ## ncRNA SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)
    
    variant.id.ncRNA <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)
    rm(variant.id.ncRNA)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))
    
    Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
    Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
    Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))
    
    rm(GENCODE.Info)
    rm(GENCODE.Info.split)
    
    variant.id.ncRNA <- seqGetData(genofile, "variant.id")
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
    is.in <- union(is.in,which(Gene_list_3==gene_name))
    
    variant.is.in <- variant.id.ncRNA[is.in]
    
  }
  
  print(variant.is.in)
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ##Genotype
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
  
  genotype <- Geno
  
  if(dim(genotype)[2] == 1){
    return(matrix(0,nrow = dim(genotype)[1],ncol = 1))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  
  print(dim(genotype))
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  
  print(RV_label)
  
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

Train_PVals_All <- read.csv(paste0("noncoding_sig.csv"))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$Chr == chr,]

## agds dir

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  ## Null Model
  obj_nullmodel_train <- get(load(paste0(trait,"_Train_Null_Model.RData")))
  system(paste0("rm ",paste0(trait,"_Train_Null_Model.RData")))
  obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
  system(paste0("rm ",paste0(trait,"_Tune_Null_Model.RData")))
  obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))
  system(paste0("rm ",paste0(trait,"_Validation_Null_Model.RData")))
  
  obj_nullmodel <- obj_nullmodel_train
  obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)
  
  ## Parameter
  QC_label <- "annotation/info/QC_label2"
  geno_missing_imputation <- "mean"
  variant_type <- "SNV"	
  
  ## Annotation_dir
  Annotation_dir <- "annotation/info/FunctionalAnnotation"
  ## Annotation channel
  Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")
  
  G_star_gene_centric_noncoding <- list()
  
  Train_PVals_All_tmp <- Train_PVals_All[Train_PVals_All$Trait == trait,]
  
  print(head(Train_PVals_All_tmp))
  
  if(nrow(Train_PVals_All_tmp) == 0){
    write.csv(NULL,file = paste0(trait,"_G_Star_Noncoding_Chr",chr,".csv"))
  }else{
    gds.path <- list.files()[str_detect(list.files(),".gds")]
    genofile <- seqOpen(gds.path)
    
    for(i in 1:nrow(Train_PVals_All)){
      ## Chr
      chr <- Train_PVals_All$Chr[i]
      ## Gene name
      gene_name <- Train_PVals_All$Gene[i]
      ## Coding mask
      category <- Train_PVals_All$Category[i]
      
      print(gene_name)
      print(category)
      
      G_star_gene_centric_noncoding[[i]] <- Gene_Centric_Noncoding_G_Star(chr=chr,gene_name=gene_name,category=category ,
                                                                          genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE)  
    }
    seqClose(genofile)
    
    G_star_gene_centric_noncoding <- do.call(cbind,G_star_gene_centric_noncoding)
    
    fwrite(G_star_gene_centric_noncoding,file = paste0(trait,"_G_Star_Noncoding_Chr",chr,".csv"),row.names = FALSE)
  } 
}

system("rm Annotation_name_catalog.csv")
system(paste0("rm ",gds.path))
system("rm coding_sig.csv")
system("rm noncoding_sig.csv")


