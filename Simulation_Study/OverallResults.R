rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- 1

results <- NULL

for(i in 1:length(Y_tune)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All",i,".RData"))
  sl_all <- SL.result
  sl_all$method <- "SL_All"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result",i,".RData"))
  ct <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2_result",i,".RData"))
  lassosum2 <- lassosum2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result",i,".RData"))
  ldpred2 <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result",i,".RData"))
  GeneCentric_Coding_STAARO <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result",i,".RData"))
  GeneCentric_Coding_Burden <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result",i,".RData"))
  GeneCentric_Noncoding_STAARO <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result",i,".RData"))
  GeneCentric_Noncoding_Burden <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result",i,".RData"))
  SlidingWindow_STAARO <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result",i,".RData"))
  SlidingWindow_Burden <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO",i,".RData"))
  BestAll_RareVariant_STAARO <- SL.result
  BestAll_RareVariant_STAARO$method <- "SL_All_STAARO"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden",i,".RData"))
  BestAll_RareVariant_Burden <- SL.result
  BestAll_RareVariant_Burden$method <- "SL_All_Burden"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result",i,".RData"))
  BestAll_STAARO <- SL.result
  BestAll_STAARO$method <- "CV_RV_STAARO_All"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result",i,".RData"))
  BestAll_Burden <- SL.result
  BestAll_Burden$method <- "CV_RV_Burden_All"
  
  
  results_tmp <- rbind(sl_all,ct,lassosum2,ldpred2,
                   GeneCentric_Coding_STAARO,GeneCentric_Coding_Burden,
                   GeneCentric_Noncoding_STAARO,GeneCentric_Noncoding_Burden,
                   SlidingWindow_STAARO,SlidingWindow_Burden,
                   BestAll_RareVariant_STAARO,BestAll_RareVariant_Burden,BestAll_Burden,BestAll_STAARO)
  
  results <- rbind(results,results_tmp)
  
  rm(list=setdiff(ls(), c("results_tmp","results","i"))) 
}

aggregate(.~method,data = results,mean)
