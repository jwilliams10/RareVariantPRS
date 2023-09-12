Burden_Effect_Size <- function(genotype,obj_nullmodel,
                  rare_maf_cutoff=0.01,rv_num_cutoff=2){

  if(class(genotype) != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

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

