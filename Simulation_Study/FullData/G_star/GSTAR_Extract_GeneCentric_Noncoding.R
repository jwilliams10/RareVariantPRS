rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricNoncoding/chr22.Rdata")
noncoding_sig <- NULL
for(i in 1:length(results_noncoding)){
  noncoding_sig <- rbind(noncoding_sig,unlist(results_noncoding[[i]][,1:3]))
}

noncoding_sig <- as.data.frame(noncoding_sig)
noncoding_sig$Chr <- as.numeric(noncoding_sig$Chr)


source("~/RareVariantPRS/Simulation_Study/FullData/G_star/g_star_gene_centric_noncoding.R")

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_agds_dir.Rdata"))

## Null Model
obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_Annotation_name_catalog.Rdata"))

obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

G_star_gene_centric_noncoding <- list()

for(i in 1:nrow(noncoding_sig)){
  chr <- noncoding_sig$Chr[i]
  ## Gene name
  gene_name <- noncoding_sig$Gene[i]
  ## Coding mask
  category <- noncoding_sig$Category[i]
  
  G_star_gene_centric_noncoding[[i]] <- Gene_Centric_Noncoding_G_Star(chr=chr,gene_name=gene_name,category=category ,
                                  genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE) 
}

G_star_gene_centric_noncoding <- do.call(cbind,G_star_gene_centric_noncoding)

seqClose(genofile)

save(G_star_gene_centric_noncoding,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Noncoding.RData")