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
library(stringr)

# for trait in 1 2 3 4 5 6;
# do
# for array in {1..1600};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/GeneCentric_Noncoding_EffectSizes.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/GeneCentric_Noncoding_EffectSizes.sh -icmd="bash GeneCentric_Noncoding_EffectSizes.sh ${trait} ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding_es/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2
# done
# done

## source code
Burden_Effect_Size <- function(genotype,obj_nullmodel,
                               rare_maf_cutoff=0.01,rv_num_cutoff=2){
  
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
    
    return(list(num_variant = num_variant,
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

Gene_Centric_Noncoding_Burden_Effect_Size_Jake <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
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
  
  burden_eff <- 0
  try(burden_eff <- Burden_Effect_Size(genotype=Geno,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  
  results <- c()
  if(class(burden_eff)=="list"){
    results_temp <- genes[1,]
    results_temp[3] <- category
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- burden_eff$num_variant
    
    results_temp <- c(results_temp,burden_eff$cMAC,burden_eff$Burden_Score_Stat,burden_eff$Burden_SE_Score,burden_eff$Burden_pvalue,burden_eff$Burden_Est,burden_eff$Burden_SE_Est)
    
    results <- rbind(results,results_temp)
  }
  
  if(!is.null(results)){
    colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
    colnames(results) <- c("Gene name","Chr","Category","#SNV",
                           "cMAC","Burden_Score_Stat","Burden_SE_Score",
                           "Burden_pvalue","Burden_Est","Burden_SE_Est")
  }
  
  seqResetFilter(genofile)
  return(results)
}

###########################################################
#           User Input
###########################################################

trait <- gsub("_noncoding_sig.csv","",list.files()[str_detect(list.files(),"_noncoding_sig.csv")])

### Significant Results 
noncoding_sig <- read_csv(paste0(trait,"_noncoding_sig.csv"))
colnames(noncoding_sig) <- c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAAR_O")

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

chunks <- list()

set <- c(0,184,288,384,448,528,616,688,744,816,880,984,1072,1096,1152,1200,1264,1360,1384,1504,1552,1568,1600)

for(i in 1:22){
  indexes <- which(noncoding_sig$Chr == i)
  a <- cut(seq_along(indexes), set[i + 1] - set[i], labels = FALSE)
  chunks <- c(chunks,split(indexes,a))
}

arrayid <- as.numeric(commandArgs(TRUE)[1])

noncoding_sig <- noncoding_sig[chunks[[arrayid]],]

chr <- unique(noncoding_sig$Chr)

gds.path <- list.files()[str_detect(list.files(),".gds")]
genofile <- seqOpen(gds.path)

unique_ids <- paste0(noncoding_sig$Gene,"__",noncoding_sig$Category)

noncoding_effectsizes_parlapply <- function(x,
                                            chr, genofile, obj_nullmodel, QC_label, variant_type, geno_missing_imputation, Annotation_dir, Annotation_name_catalog){
  tmp <- unlist(strsplit(x,"__"))
  gene <- tmp[1]
  category <- tmp[2]
  
  a <- Gene_Centric_Noncoding_Burden_Effect_Size_Jake(chr = chr,gene_name = gene,category=category,
    genofile = genofile,obj_nullmodel = obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
    Annotation_dir=Annotation_dir,Annotation_name_catalog = Annotation_name_catalog,silent=FALSE)
  
  if(!is.null(a)){
    a <- data.frame(Gene = gene,Chr = chr,Category = category,Number_SNV = a[[4]],cMAC = a[[5]],Burden_Score_Stat = a[[6]],Burden_SE_Score = a[[7]],Burden_pvalue = a[[8]],Burden_Est = a[[9]], Burden_SE_Est = a[[10]])
    a <- unique(a)
  }
  return(a)
} 

a <- lapply(unique_ids,noncoding_effectsizes_parlapply,
            chr = chr, genofile = genofile, obj_nullmodel = obj_nullmodel, QC_label = QC_label, variant_type = variant_type, geno_missing_imputation = geno_missing_imputation, Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog)

seqClose(genofile)

effect_sizes <- data.table::rbindlist(a)
effect_sizes$Chr <- unique(noncoding_sig$Chr)

noncoding_sig <- inner_join(noncoding_sig,effect_sizes)

write.csv(noncoding_sig,row.names = FALSE,file = paste0(trait,"_noncoding_sig_array",arrayid,".csv"))

a <- list.files()[list.files() != paste0(trait,"_noncoding_sig_array",arrayid,".csv")]

for(i in 1:length(a)){
  system(paste0("rm ",a[i]))
}
