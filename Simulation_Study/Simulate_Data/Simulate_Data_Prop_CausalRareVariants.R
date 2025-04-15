rm(list = ls())
set.seed(1340)
gc()

library(bigsnpr)
library(dplyr)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(STAAR)
library(STAARpipeline)

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

if(file.exists("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")){
  system("rm /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")
  system("rm /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.bk")
}

rds <- bigsnpr::snp_readBed("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.bed") 

# Loading the data from backing files
common_variants <- snp_attach("/gpfs/gsfs12/users/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")

fam_file <- read.table("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam", quote="\"", comment.char="")

number_snps <- dim(common_variants$genotypes)[2]

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Coding.RData")

obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

id.genotype <- seqGetData(genofile,"sample.id")
id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
phenotype.id <- as.character(obj_nullmodel$id_include)
phenotype.id.merge <- data.frame(phenotype.id)
phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))

seqClose(genofile)

number_genes <- dim(G_star_gene_centric_coding)[2]

causalprop_vec <- c(0.2,0.05,0.01,0.001,0.0005)
scale <- c(0,1)

Y <- list()

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

h2_overall_EUR <- vector()
h2_overall_EAS <- vector()
h2_overall_SAS <- vector()
h2_overall_AMR <- vector()
h2_overall_AFR <- vector()

h2_common_EUR <- vector()
h2_common_EAS <- vector()
h2_common_SAS <- vector()
h2_common_AMR <- vector()
h2_common_AFR <- vector()

MAF_common_EUR <- vector()
MAF_common_EAS <- vector()
MAF_common_SAS <- vector()
MAF_common_AMR <- vector()
MAF_common_AFR <- vector()

Var_common_EUR <- vector()
Var_common_EAS <- vector()
Var_common_SAS <- vector()
Var_common_AMR <- vector()
Var_common_AFR <- vector()

Average_Burden_rare_EUR <- vector()
Average_Burden_rare_EAS <- vector()
Average_Burden_rare_SAS <- vector()
Average_Burden_rare_AMR <- vector()
Average_Burden_rare_AFR <- vector()

Var_Burden_rare_EUR <- vector()
Var_Burden_rare_EAS <- vector()
Var_Burden_rare_SAS <- vector()
Var_Burden_rare_AMR <- vector()
Var_Burden_rare_AFR <- vector()

h2_rare_EUR <- vector()
h2_rare_EAS <- vector()
h2_rare_SAS <- vector()
h2_rare_AMR <- vector()
h2_rare_AFR <- vector()

count <- 1

h2_common <- 0.05
h2_rare <- 0.05/4

for(j in 1:length(causalprop_vec)){
  number_causal_snps <- round(causalprop_vec[j]*number_snps)
  causal_snps <- sample(1:number_snps,number_causal_snps)
  
  beta_snps <- rnorm(number_causal_snps,mean = 0,sqrt(h2_common/number_causal_snps))
  
  g_snps <- common_variants$genotypes[,causal_snps]
  for(i in 1:ncol(g_snps)){
    g_snps[is.na(g_snps[,i]),i] <- 0
  }
  scaled_causal_snps <- scale(g_snps)
  
  number_causal_genes <- round(causalprop_vec[j]*number_genes)
  causal_genes <- sample(1:number_genes,number_causal_genes)
  
  load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricCoding/chr22.Rdata")
  coding_sig <- NULL
  for(i in 1:length(results_coding)){
    coding_sig <- rbind(coding_sig,unlist(results_coding[[i]][,1:3]))
  }
  coding_sig <- as.data.frame(coding_sig)
  coding_sig$Chr <- as.numeric(coding_sig$Chr)
  coding_sig <- coding_sig[causal_genes,]
  
  QC_label <- "annotation/info/QC_label"
  geno_missing_imputation <- "mean"
  variant_type <- "SNV"	
  
  ## agds dir
  agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))
  ## Null Model
  obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))
  ## Annotation_dir
  Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
  ## Annotation channel
  Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_Annotation_name_catalog.Rdata"))
  
  genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")
  
  G_star_gene_centric_coding <- list()
  for(rows in 1:nrow(coding_sig)){
    chr <- coding_sig$Chr[rows]
    ## Gene name
    gene_name <- coding_sig$Gene[rows]
    ## Coding mask
    category <- coding_sig$Category[rows]
    
    rarevariants_gene_centric_coding <- Gene_Centric_Coding_RareVariants(chr=chr,gene_name=gene_name,category=category ,
                                                                              genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE) 
    rarevariants_probability_causal <- runif(1,min = 0.1,max = 0.9)
    number_causal_rarevariants <- ceiling(ncol(rarevariants_gene_centric_coding)*rarevariants_probability_causal)
    causal_rarevariants <- sample(1:ncol(rarevariants_gene_centric_coding),number_causal_rarevariants)
    G_star_gene_centric_coding[[rows]] <- rarevariants_gene_centric_coding[,causal_rarevariants]%*%matrix(1,nrow=number_causal_rarevariants,ncol = 1)
  }
  
  seqClose(genofile)
  
  G_star_gene_centric_coding <- do.call(cbind,G_star_gene_centric_coding)
  
  scaled_causal_genes <- scale(G_star_gene_centric_coding)
  
  beta_genes <- rnorm(number_causal_genes,mean = 0,sqrt(h2_rare/number_causal_genes))
  
  for(q in 1:length(scale)){
    for(l in 1:100){
      
      
      if(scale[q] == 0){
        # Var(G Beta) = c c/v -> 0.05
        v <- var(g_snps%*%matrix(beta_snps,ncol = 1))/h2_common
        Y_hat_1 <- data.frame(IDs = fam_file[,2],Y_hat_Common = g_snps%*%matrix(beta_snps/sqrt(as.numeric(v)),ncol = 1))
        
      }else{
        v <- var(scaled_causal_snps%*%matrix(beta_snps,ncol = 1))/h2_common
        Y_hat_1 <- data.frame(IDs = fam_file[,2],Y_hat_Common = scaled_causal_snps%*%matrix(beta_snps/sqrt(as.numeric(v)),ncol = 1)) 
      }
      
      if(scale[q] == 0){
        # Var(G Beta) = c c/v -> 0.05
        v <- var(G_star_gene_centric_coding%*%matrix(beta_genes,ncol = 1))/h2_rare
        Y_hat_2 <- data.frame(IDs = as.numeric(phenotype.id.merge[,1]),Y_hat_Rare = G_star_gene_centric_coding%*%matrix(beta_genes/sqrt(as.numeric(v)),ncol = 1))
      }else{
        v <- var(scaled_causal_genes%*%matrix(beta_genes,ncol = 1))/h2_rare
        Y_hat_2 <- data.frame(IDs = as.numeric(phenotype.id.merge[,1]),Y_hat_Rare = scaled_causal_genes%*%matrix(beta_genes/sqrt(as.numeric(v)),ncol = 1)) 
      }
      
      ### Combine
      
      Y_hat <- inner_join(Y_hat_1,Y_hat_2)
      
      epsilon <- rnorm(nrow(Y_hat),mean = 0,sd = sqrt(1 - h2_common - h2_rare))
      
      Y_raw <- Y_hat$Y_hat_Common + Y_hat$Y_hat_Rare
      
      Y_tmp <- Y_raw + epsilon
      
      h2_common_EUR[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])
      h2_common_SAS[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])
      h2_common_AFR[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])
      h2_common_AMR[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])
      h2_common_EAS[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])
      
      h2_rare_EUR[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])
      h2_rare_SAS[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])
      h2_rare_AFR[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])
      h2_rare_AMR[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])
      h2_rare_EAS[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])
      
      MAF_common_EUR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_SAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_AFR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_AMR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_EAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],],2,function(x){sum(x)/(2*length(x))}))
      
      Var_common_EUR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],],2,var))
      Var_common_SAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],],2,var))
      Var_common_AFR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],],2,var))
      Var_common_AMR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],],2,var))
      Var_common_EAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],],2,var))
      
      Average_Burden_rare_EUR[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],,drop = FALSE]))
      Average_Burden_rare_SAS[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],,drop = FALSE]))
      Average_Burden_rare_AFR[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],,drop = FALSE]))
      Average_Burden_rare_AMR[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],,drop = FALSE]))
      Average_Burden_rare_EAS[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],,drop = FALSE]))
      
      Var_Burden_rare_EUR[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],,drop = FALSE],2,var))
      Var_Burden_rare_SAS[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],,drop = FALSE],2,var))
      Var_Burden_rare_AFR[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],,drop = FALSE],2,var))
      Var_Burden_rare_AMR[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],,drop = FALSE],2,var))
      Var_Burden_rare_EAS[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],,drop = FALSE],2,var))
      
      h2_overall_EUR[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])
      h2_overall_SAS[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])
      h2_overall_AFR[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])
      h2_overall_AMR[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])
      h2_overall_EAS[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])
      
      Y[[count]] <- data.frame(IDs = Y_hat$IDs,Y=Y_tmp)
      
      print(count)
      
      count <- count + 1
    }
  }
  gc()
}

h2_dat <- data.frame(causal_prop = rep(rep(causalprop_vec,each = 200),5),scaled = rep(rep(rep(c("No","Yes"),each = 100),5),5),
                     Ancestry = rep(c("EUR","SAS","AFR","AMR","EAS"),each = 1000),
                     h2_common = c(h2_common_EUR,h2_common_SAS,h2_common_AFR,h2_common_AMR,h2_common_EAS),
                     h2_rare = c(h2_rare_EUR,h2_rare_SAS,h2_rare_AFR,h2_rare_AMR,h2_rare_EAS),
                     h2_overall = c(h2_overall_EUR,h2_overall_SAS,h2_overall_AFR,h2_overall_AMR,h2_overall_EAS),
                     MAF_common = c(MAF_common_EUR,MAF_common_SAS,MAF_common_AFR,MAF_common_AMR,MAF_common_EAS),
                     Var_common = c(Var_common_EUR,Var_common_SAS,Var_common_AFR,Var_common_AMR,Var_common_EAS),
                     Average_Burden_rare = c(Average_Burden_rare_EUR,Average_Burden_rare_SAS,Average_Burden_rare_AFR,Average_Burden_rare_AMR,Average_Burden_rare_EAS),
                     Var_Burden_rare = c(Var_Burden_rare_EUR,Var_Burden_rare_SAS,Var_Burden_rare_AFR,Var_Burden_rare_AMR,Var_Burden_rare_EAS))

write.csv(h2_dat,file = "/data/williamsjacr/UKB_WES_Simulation/Sim_Characteristics_PropRareVariants.csv",row.names = FALSE)

lapply(Y,function(x){var(x$Y)})
summary(unlist(lapply(Y,function(x){var(x$Y)})))

save(Y,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/Y_n_",nrow(Y_hat),"_h2_common_",h2_common,"_h2_rare_",h2_rare,".RData"))
save(Y,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/Y_n_",nrow(Y_hat),"_h2_common_",h2_common,"_h2_rare_",h2_rare,".RData"))
