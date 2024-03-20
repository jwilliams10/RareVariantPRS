rm(list=ls())
list.files()
gc()

# for trait in 1 2 3 4 5;
# do
# for array in {1..1172};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/GeneCentric_ncRNA_PRS_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/GeneCentric_ncRNA_PRS_Binary.sh -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/ncRNA_chr.csv -icmd="bash GeneCentric_ncRNA_PRS_Binary.sh ${trait} ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentric_ncRNA_PRS/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2
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

Gene_Centric_Noncoding_Burden_PRS <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
                                              genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                              BETA,
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
  
  Burden_PRS <- Burden_PRS(genotype=Geno,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,BETA = BETA)
  
  Burden_PRS <- data.frame(ID = phenotype.id.merge$phenotype.id,PRS = Burden_PRS)
  
  seqResetFilter(genofile)
  return(Burden_PRS)
}


###########################################################
#           User Input
###########################################################

trait <- gsub("_ncRNA_sig_es.csv","",list.files()[str_detect(list.files(),"_ncRNA_sig_es.csv")])

### Significant Results 
Train_Effect_Sizes_All <- read_csv(paste0(trait,"_ncRNA_sig_es.csv"))

## Null model
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

if(arrayid>586){
  Burden <- 1
  arrayid <- arrayid - 586
}else{
  Burden <- 0
}

sets <- c(0,22,44,66,88,110,132,154,176,198,220,242,264,286,386,586)
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
}else if(arrayid <= 286){
  threshold <- 13
  jobid <- arrayid - 264
}else if(arrayid <= 386){
  threshold <- 14
  jobid <- arrayid - 286
  set <- c(0,9,20,27,32,37,42,46,50,54,58,62,66,71,74,78,82,86,89,93,96,98,100)
}else{
  threshold <- 15
  jobid <- arrayid - 386
  set <- c(0,18,40,55,64,75,84,93,101,109,116,125,133,142,149,157,164,173,178,186,192,197,200)
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
      a <- cut(seq_along(indexes), set[i + 1] - set[i], labels = FALSE)
      chunks <- c(chunks,split(indexes,a))
    }
    Train_Effect_Sizes_All <- Train_Effect_Sizes_All[chunks[[jobid]],]
  }else{
    Train_Effect_Sizes_All <- Train_Effect_Sizes_All[Train_Effect_Sizes_All$Chr == jobid,]
  }
  
  if(nrow(Train_Effect_Sizes_All) == 0){
    PRS <- data.frame(ID = 1:length(obj_nullmodel$id_include),PRS = 0)
  }else{
    chr <- unique(Train_Effect_Sizes_All$Chr)
    
    gds.path <- list.files()[str_detect(list.files(),".gds")]
    genofile <- seqOpen(gds.path)
    
    unique_ids <- paste0(Train_Effect_Sizes_All$Gene,"__",Train_Effect_Sizes_All$Category,"__",Train_Effect_Sizes_All$Burden_Est)
    
    ncRNA_effectsizes_parlapply <- function(x,
                                            chr, genofile, obj_nullmodel, QC_label, variant_type, geno_missing_imputation, Annotation_dir, Annotation_name_catalog){
      
      tmp <- unlist(strsplit(x,"__"))
      gene <- tmp[1]
      category <- tmp[2]
      BETA <- as.numeric(tmp[3])
      
      a <- Gene_Centric_Noncoding_Burden_PRS(chr = chr,gene_name = gene,category=category, BETA = BETA,
                                             genofile = genofile,obj_nullmodel = obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog = Annotation_name_catalog,silent=FALSE)
      return(a)
    } 
    
    PRSs <- lapply(unique_ids,ncRNA_effectsizes_parlapply,
                   chr = chr, genofile = genofile, obj_nullmodel = obj_nullmodel, QC_label = QC_label, variant_type = variant_type, geno_missing_imputation = geno_missing_imputation, Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog)
    
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

write.csv(PRS,row.names = FALSE,file = paste0(trait,"_ncRNA_PRS",arrayid_original,".csv"))

a <- list.files()[list.files() != paste0(trait,"_ncRNA_PRS",arrayid_original,".csv")]

for(i in 1:length(a)){
  system(paste0("rm ",a[i]))
}
