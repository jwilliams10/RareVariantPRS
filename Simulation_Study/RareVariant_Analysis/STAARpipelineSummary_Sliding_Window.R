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
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_Analysis",i,".Rdata"))
  
  slidingwindow_sig <- data.frame(Chr = unlist(results_sliding_window[,1]),Start = unlist(results_sliding_window[,2]),End = unlist(results_sliding_window[,3]),Burden_1_1 = unlist(results_sliding_window[,47]),STAARO = unlist(results_sliding_window[,90]))
  
  slidingwindow_sig$Chr <- as.numeric(slidingwindow_sig$Chr) 
  
  save(slidingwindow_sig,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/slidingwindow_sig",i,".Rdata"))
}