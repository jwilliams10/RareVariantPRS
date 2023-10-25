##########################################################
# Summarization and visualization of gene-centric 
# noncoding analysis results using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
##########################################################
rm(list=ls())
gc()

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

for(i in 1:length(Y_tune)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/Train_Analysis",i,".Rdata"))
  noncoding_sig <- NULL
  for(j in 1:length(results_noncoding)){
    noncoding_sig <- rbind(noncoding_sig,unlist(results_noncoding[[j]][,c(1:3,47,90)]))
  }
  
  noncoding_sig <- as.data.frame(noncoding_sig)
  noncoding_sig$Chr <- as.numeric(noncoding_sig$Chr) 
  
  save(noncoding_sig,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/noncoding_sig",i,".Rdata"))
}