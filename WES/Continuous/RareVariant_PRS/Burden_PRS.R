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

