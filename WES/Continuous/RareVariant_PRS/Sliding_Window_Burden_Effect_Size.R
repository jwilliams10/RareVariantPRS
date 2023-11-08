Sliding_Window_Burden_Effect_Size <- function(chr,genofile,obj_nullmodel,start_loc,end_loc,
											rare_maf_cutoff=0.01,rv_num_cutoff=2,
											QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
											silent=FALSE){

	## evaluate choices
	geno_missing_imputation <- match.arg(geno_missing_imputation)
	variant_type <- match.arg(variant_type)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")

	## Position
	position <- as.numeric(seqGetData(genofile, "position"))

	is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
    seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")

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

	burden_eff <- 0
	try(burden_eff <- Burden_Effect_Size(genotype=Geno,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results <- c()
	if(class(burden_eff)=="list"){
		results_temp <- c(chr,start_loc,end_loc,burden_eff$num_variant)
		results_temp <- c(results_temp,burden_eff$cMAC,burden_eff$Burden_Score_Stat,burden_eff$Burden_SE_Score,burden_eff$Burden_pvalue,burden_eff$Burden_Est,burden_eff$Burden_SE_Est)

		results <- rbind(results,results_temp)
	}

	if(!is.null(results)){
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results) <- c("Chr","Start Loc","End Loc","#SNV",
								"cMAC","Burden_Score_Stat","Burden_SE_Score",
								"Burden_pvalue","Burden_Est","Burden_SE_Est")
	}
	seqResetFilter(genofile)
	
	return(results)
}
