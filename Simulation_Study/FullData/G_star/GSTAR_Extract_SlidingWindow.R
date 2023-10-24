rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/SlidingWindow/chr22.Rdata")
sliding_window_sig <- results_sliding_window[,1:3]

sliding_window_sig <- as.data.frame(sliding_window_sig)
sliding_window_sig$Chr <- as.numeric(sliding_window_sig$Chr)

sliding_window_sig$Start_Loc <- unlist(sliding_window_sig$`Start Loc`)
sliding_window_sig$End_Loc <- unlist(sliding_window_sig$`End Loc`)

sliding_window_sig <- sliding_window_sig[,c(1,4,5)]


source("~/RareVariantPRS/Simulation_Study/FullData/G_star/g_star_sliding_window.R")

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

G_star_sliding_window <- list()

for(i in 1:nrow(sliding_window_sig)){
  chr <- sliding_window_sig$Chr[i]
  ## Gene name
  start <- sliding_window_sig$Start_Loc[i]
  ## Coding mask
  end <- sliding_window_sig$End_Loc[i]
  
  G_star_sliding_window[[i]] <- SlidingWindow_G_Star(chr=chr,start_loc = start, end = end ,
                                     genofile,obj_nullmodel = obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,silent=FALSE) 
}

G_star_sliding_window <- do.call(cbind,G_star_sliding_window)

seqClose(genofile)

save(G_star_sliding_window,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/SlidingWindow.RData")