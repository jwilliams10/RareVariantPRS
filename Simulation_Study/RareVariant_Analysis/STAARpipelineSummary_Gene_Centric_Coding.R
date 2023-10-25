##########################################################
# Summarization and visualization of gene-centric 
# coding analysis results using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
##########################################################
rm(list=ls())
gc()

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

for(i in 1:length(Y_tune)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/Train_Analysis",i,".Rdata"))
  coding_sig <- NULL
  for(j in 1:length(results_coding)){
    coding_sig <- rbind(coding_sig,unlist(results_coding[[j]][,c(1:3,47,90)]))
  }
  
  coding_sig <- as.data.frame(coding_sig)
  coding_sig$Chr <- as.numeric(coding_sig$Chr) 
  
  save(coding_sig,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/coding_sig",i,".Rdata"))
}

