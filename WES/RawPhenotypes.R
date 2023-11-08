rm(list = ls())

library(data.table)
library(dplyr)

ukb_pheno <- fread("/data/BB_Bioinformatics/ProjectData/UKB/Processed_phenotype_data/ukb_multi.pheno")
ukb_covariate <- fread("/data/BB_Bioinformatics/ProjectData/UKB/Processed_phenotype_data/ukb_multi.cov")

ukb_covariate <- ukb_covariate[,c("FID","IID","ancestry","female","age",paste0("pc",1:10))]
ukb_covariate <- ukb_covariate[ukb_covariate$ancestry == "EUR",]
colnames(ukb_covariate) <- c(colnames(ukb_covariate)[1:3],"sex",colnames(ukb_covariate)[5:15])
ukb_covariate <- ukb_covariate[,-c(3)]

ukb_pheno <- inner_join(ukb_pheno,ukb_covariate)

rm(ukb_covariate)

save(ukb_pheno,file = "/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
