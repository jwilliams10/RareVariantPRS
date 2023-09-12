Gene_Centric_Noncoding_Burden_Effect_Size <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
                                   genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
								   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){

	## evaluate choices
	category <- match.arg(category)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)
	
	genes <- genes_info[genes_info[,2]==chr,]

	if(category=="downstream")
	{
		results <- downstream_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	if(category=="upstream")
	{
		results <- upstream_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	if(category=="UTR")
	{
		results <- UTR_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	if(category=="promoter_CAGE")
	{
		results <- promoter_CAGE_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	if(category=="promoter_DHS")
	{
		results <- promoter_DHS_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	if(category=="enhancer_CAGE")
	{
		results <- enhancer_CAGE_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	if(category=="enhancer_DHS")
	{
		results <- enhancer_DHS_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}
	
	if(category=="ncRNA")
	{
		results <- ncRNA_burden_effect_size(chr,gene_name,genofile,obj_nullmodel,genes,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=silent)
	}

	return(results)
}

