rm(list = ls())

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/GeneCentricNoncoding/chr22.Rdata")
noncoding_sig <- NULL
for(i in 1:length(results_noncoding)){
  noncoding_sig <- rbind(noncoding_sig,unlist(results_noncoding[[i]][,1:3]))
}

noncoding_sig <- as.data.frame(noncoding_sig)
noncoding_sig$Chr <- as.numeric(noncoding_sig$Chr)
noncoding_sig_overall <- noncoding_sig
colnames(noncoding_sig_overall)[1] <- "Gene"

noncoding_sig <- read.csv("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes1.csv")
noncoding_sig <- left_join(noncoding_sig_overall,noncoding_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Tune_Null_Model1.RData"))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Noncoding.RData")

G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[,!is.na(noncoding_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

i <- as.numeric(commandArgs(TRUE)[1])

noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_EffectSizes",i,".csv"))
prs_mat <- NULL

for(j in 1:2){
  if(j == 1){
    p_values <- noncoding_sig$Burden_pvalue
  }else{
    p_values <- noncoding_sig$STAAR_O
  }
  for(k in 1:length(thresholds)){
    beta <- matrix(0,ncol = 1,nrow = ncol(G_star_gene_centric_noncoding))
    beta[p_values < thresholds[k]] <- noncoding_sig$Burden_Est[p_values < thresholds[k]]
    prs_mat <- cbind(prs_mat,G_star_gene_centric_noncoding%*%beta)
  }
}
save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/tune_prs_mat",i,".RData"))



