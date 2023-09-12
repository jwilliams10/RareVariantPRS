promoter_DHS_burden_effect_size <- function(chr,gene_name,genofile,obj_nullmodel,genes,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

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
	gc()

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
	variant.id.SNV <- variant.id[SNVlist]

	dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
	dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
	dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
	dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

	seqResetFilter(genofile)

	rm(dfPromrOCRsVarGene)
	gc()

	### Gene
	is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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

	burden_eff <- 0
	try(burden_eff <- Burden_Effect_Size(genotype=Geno,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results <- c()
	if(class(burden_eff)=="list")
	{
		results_temp <- genes[1,]
		results_temp[3] <- "promoter_DHS"
		results_temp[2] <- chr
		results_temp[1] <- as.character(gene_name)
		results_temp[4] <- burden_eff$num_variant

		results_temp <- c(results_temp,burden_eff$cMAC,burden_eff$Burden_Score_Stat,burden_eff$Burden_SE_Score,burden_eff$Burden_pvalue,burden_eff$Burden_Est,burden_eff$Burden_SE_Est)

		results <- rbind(results,results_temp)
	}

	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results) <- c("Gene name","Chr","Category","#SNV",
								"cMAC","Burden_Score_Stat","Burden_SE_Score",
								"Burden_pvalue","Burden_Est","Burden_SE_Est")
	}

	seqResetFilter(genofile)
	return(results)
}

