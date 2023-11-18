rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")

results_70_EUR <- NULL
results_70_NonEur <- NULL
results_70_EAS <- NULL
results_70_AFR <- NULL
results_70_SAS <- NULL
results_70_MIX <- NULL
results_70_UNK <- NULL

for(i in 1:length(Y_train)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_Eur",i,".RData"))
  sl_all_Eur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_NonEur",i,".RData"))
  sl_all_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_EAS",i,".RData"))
  sl_all_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_AFR",i,".RData"))
  sl_all_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_SAS",i,".RData"))
  sl_all_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_MIX",i,".RData"))
  sl_all_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All_UNK",i,".RData"))
  sl_all_UNK <- SL.result
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_EUR",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_EUR <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_NonEur",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_NonEur <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_EAS",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_EAS <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_AFR",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_AFR <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_SAS",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_SAS <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_MIX",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_MIX <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result_UNK",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_UNK <- CT.result[[1]]
  
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_EUR",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_EUR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_NonEur",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_NonEur <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_EAS",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_EAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_AFR",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_AFR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_SAS",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_SAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_MIX",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_MIX <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/LASSOSUM2_result_UNK",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_UNK <- LASSOSUM2.result[[1]]
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_EUR",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_EUR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_NonEur",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_NonEur <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_EAS",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_EAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_AFR",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_AFR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_SAS",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_SAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_MIX",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_MIX <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result_UNK",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_UNK <- ldpred2.result[[1]]
  
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_EUR",i,".RData"))
  GeneCentric_Coding_STAARO_EUR <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_NonEUR",i,".RData"))
  GeneCentric_Coding_STAARO_NonEur <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_EAS",i,".RData"))
  GeneCentric_Coding_STAARO_EAS <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_AFR",i,".RData"))
  GeneCentric_Coding_STAARO_AFR <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_SAS",i,".RData"))
  GeneCentric_Coding_STAARO_SAS <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_MIX",i,".RData"))
  GeneCentric_Coding_STAARO_MIX <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_UNK",i,".RData"))
  GeneCentric_Coding_STAARO_UNK <- r2.result_GeneCentric_Coding_STAARO
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_EUR",i,".RData"))
  GeneCentric_Coding_Burden_EUR <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_NonEUR",i,".RData"))
  GeneCentric_Coding_Burden_NonEur <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_EAS",i,".RData"))
  GeneCentric_Coding_Burden_EAS <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_AFR",i,".RData"))
  GeneCentric_Coding_Burden_AFR <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_SAS",i,".RData"))
  GeneCentric_Coding_Burden_SAS <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_MIX",i,".RData"))
  GeneCentric_Coding_Burden_MIX <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_UNK",i,".RData"))
  GeneCentric_Coding_Burden_UNK <- r2.result_GeneCentric_Coding_Burden
  
  
  

  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_EUR",i,".RData"))
  GeneCentric_Noncoding_STAARO_EUR <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_NonEUR",i,".RData"))
  GeneCentric_Noncoding_STAARO_NonEur <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_EAS",i,".RData"))
  GeneCentric_Noncoding_STAARO_EAS <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_AFR",i,".RData"))
  GeneCentric_Noncoding_STAARO_AFR <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_SAS",i,".RData"))
  GeneCentric_Noncoding_STAARO_SAS <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_MIX",i,".RData"))
  GeneCentric_Noncoding_STAARO_MIX <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_UNK",i,".RData"))
  GeneCentric_Noncoding_STAARO_UNK <- r2.result_GeneCentric_Noncoding_STAARO
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_EUR",i,".RData"))
  GeneCentric_Noncoding_Burden_EUR <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_NonEUR",i,".RData"))
  GeneCentric_Noncoding_Burden_NonEur <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_EAS",i,".RData"))
  GeneCentric_Noncoding_Burden_EAS <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_AFR",i,".RData"))
  GeneCentric_Noncoding_Burden_AFR <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_SAS",i,".RData"))
  GeneCentric_Noncoding_Burden_SAS <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_MIX",i,".RData"))
  GeneCentric_Noncoding_Burden_MIX <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_UNK",i,".RData"))
  GeneCentric_Noncoding_Burden_UNK <- r2.result_GeneCentric_Noncoding_Burden
  
  
  

  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_EUR",i,".RData"))
  SlidingWindow_STAARO_EUR <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_NonEUR",i,".RData"))
  SlidingWindow_STAARO_NonEur <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_EAS",i,".RData"))
  SlidingWindow_STAARO_EAS <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_AFR",i,".RData"))
  SlidingWindow_STAARO_AFR <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_SAS",i,".RData"))
  SlidingWindow_STAARO_SAS <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_MIX",i,".RData"))
  SlidingWindow_STAARO_MIX <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result_UNK",i,".RData"))
  SlidingWindow_STAARO_UNK <- r2.result_SlidingWindow_STAARO
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_EUR",i,".RData"))
  SlidingWindow_Burden_EUR <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_NonEUR",i,".RData"))
  SlidingWindow_Burden_NonEur <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_EAS",i,".RData"))
  SlidingWindow_Burden_EAS <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_AFR",i,".RData"))
  SlidingWindow_Burden_AFR <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_SAS",i,".RData"))
  SlidingWindow_Burden_SAS <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_MIX",i,".RData"))
  SlidingWindow_Burden_MIX <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result_UNK",i,".RData"))
  SlidingWindow_Burden_UNK <- r2.result_SlidingWindow_Burden
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_EUR",i,".RData"))
  BestAll_RareVariant_STAARO_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_NonEUR",i,".RData"))
  BestAll_RareVariant_STAARO_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_EAS",i,".RData"))
  BestAll_RareVariant_STAARO_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_AFR",i,".RData"))
  BestAll_RareVariant_STAARO_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_SAS",i,".RData"))
  BestAll_RareVariant_STAARO_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_MIX",i,".RData"))
  BestAll_RareVariant_STAARO_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_UNK",i,".RData"))
  BestAll_RareVariant_STAARO_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_EUR",i,".RData"))
  BestAll_RareVariant_Burden_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_NonEUR",i,".RData"))
  BestAll_RareVariant_Burden_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_EAS",i,".RData"))
  BestAll_RareVariant_Burden_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_AFR",i,".RData"))
  BestAll_RareVariant_Burden_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_SAS",i,".RData"))
  BestAll_RareVariant_Burden_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_MIX",i,".RData"))
  BestAll_RareVariant_Burden_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden_UNK",i,".RData"))
  BestAll_RareVariant_Burden_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_EUR",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_NonEur",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_EAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_AFR",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_SAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_MIX",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result_UNK",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_UNK <- SL.result
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_EUR",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_NonEur",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_EAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_AFR",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_SAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_MIX",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result_UNK",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_UNK <- SL.result
  
  
  results_tmp_EUR <- rbind(sl_all_Eur,ct_EUR,lassosum2_EUR,ldpred2_EUR,
                       GeneCentric_Coding_STAARO_EUR,GeneCentric_Coding_Burden_EUR,
                       GeneCentric_Noncoding_STAARO_EUR,GeneCentric_Noncoding_Burden_EUR,
                       SlidingWindow_STAARO_EUR,SlidingWindow_Burden_EUR,
                       BestAll_RareVariant_STAARO_EUR,BestAll_RareVariant_Burden_EUR,BestAll_Burden_EUR,BestAll_STAARO_EUR)
  
  results_70_EUR <- rbind(results_70_EUR,results_tmp_EUR)
  
  
  results_tmp_NonEur <- rbind(sl_all_NonEur,ct_NonEur,lassosum2_NonEur,ldpred2_NonEur,
                           GeneCentric_Coding_STAARO_NonEur,GeneCentric_Coding_Burden_NonEur,
                           GeneCentric_Noncoding_STAARO_NonEur,GeneCentric_Noncoding_Burden_NonEur,
                           SlidingWindow_STAARO_NonEur,SlidingWindow_Burden_NonEur,
                           BestAll_RareVariant_STAARO_NonEur,BestAll_RareVariant_Burden_NonEur,BestAll_Burden_NonEur,BestAll_STAARO_NonEur)
  
  results_70_NonEur <- rbind(results_70_NonEur,results_tmp_NonEur)
  
  
  results_tmp_AFR <- rbind(sl_all_AFR,ct_AFR,lassosum2_AFR,ldpred2_AFR,
                              GeneCentric_Coding_STAARO_AFR,GeneCentric_Coding_Burden_AFR,
                              GeneCentric_Noncoding_STAARO_AFR,GeneCentric_Noncoding_Burden_AFR,
                              SlidingWindow_STAARO_AFR,SlidingWindow_Burden_AFR,
                              BestAll_RareVariant_STAARO_AFR,BestAll_RareVariant_Burden_AFR,BestAll_Burden_AFR,BestAll_STAARO_AFR)
  
  results_70_AFR <- rbind(results_70_AFR,results_tmp_AFR)
  
  results_tmp_SAS <- rbind(sl_all_SAS,ct_SAS,lassosum2_SAS,ldpred2_SAS,
                           GeneCentric_Coding_STAARO_SAS,GeneCentric_Coding_Burden_SAS,
                           GeneCentric_Noncoding_STAARO_SAS,GeneCentric_Noncoding_Burden_SAS,
                           SlidingWindow_STAARO_SAS,SlidingWindow_Burden_SAS,
                           BestAll_RareVariant_STAARO_SAS,BestAll_RareVariant_Burden_SAS,BestAll_Burden_SAS,BestAll_STAARO_SAS)
  
  results_70_SAS <- rbind(results_70_SAS,results_tmp_SAS)
  
  results_tmp_EAS <- rbind(sl_all_EAS,ct_EAS,lassosum2_EAS,ldpred2_EAS,
                           GeneCentric_Coding_STAARO_EAS,GeneCentric_Coding_Burden_EAS,
                           GeneCentric_Noncoding_STAARO_EAS,GeneCentric_Noncoding_Burden_EAS,
                           SlidingWindow_STAARO_EAS,SlidingWindow_Burden_EAS,
                           BestAll_RareVariant_STAARO_EAS,BestAll_RareVariant_Burden_EAS,BestAll_Burden_EAS,BestAll_STAARO_EAS)
  
  results_70_EAS <- rbind(results_70_EAS,results_tmp_EAS)
  
  results_tmp_MIX <- rbind(sl_all_MIX,ct_MIX,lassosum2_MIX,ldpred2_MIX,
                           GeneCentric_Coding_STAARO_MIX,GeneCentric_Coding_Burden_MIX,
                           GeneCentric_Noncoding_STAARO_MIX,GeneCentric_Noncoding_Burden_MIX,
                           SlidingWindow_STAARO_MIX,SlidingWindow_Burden_MIX,
                           BestAll_RareVariant_STAARO_MIX,BestAll_RareVariant_Burden_MIX,BestAll_Burden_MIX,BestAll_STAARO_MIX)
  
  results_70_MIX <- rbind(results_70_MIX,results_tmp_MIX)
  
  results_tmp_UNK <- rbind(sl_all_UNK,ct_UNK,lassosum2_UNK,ldpred2_UNK,
                           GeneCentric_Coding_STAARO_UNK,GeneCentric_Coding_Burden_UNK,
                           GeneCentric_Noncoding_STAARO_UNK,GeneCentric_Noncoding_Burden_UNK,
                           SlidingWindow_STAARO_UNK,SlidingWindow_Burden_UNK,
                           BestAll_RareVariant_STAARO_UNK,BestAll_RareVariant_Burden_UNK,BestAll_Burden_UNK,BestAll_STAARO_UNK)
  
  results_70_UNK <- rbind(results_70_UNK,results_tmp_UNK)
  
  
  rm(list=setdiff(ls(), c("results_70_EUR","results_70_NonEur","results_70_EAS","results_70_AFR","results_70_SAS","results_70_MIX","results_70_UNK",
                          "i","Y_train"))) 
}

scale <- rep(rep(c("Unscaled","Scaled"),each = 20*14),5)
causal_prop <- rep(c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01","Causal Prop. 0.001","Causal Prop. 0.0005"),each = 40*14)

results_70_EUR <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_EUR$method,R2 = results_70_EUR$r2,R2_Low = results_70_EUR$r2_low,R2_High = results_70_EUR$r2_high)
results_70_EUR$Train_Size <- nrow(Y_train[[1]])

results_70_NonEur <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_NonEur$method,R2 = results_70_NonEur$r2,R2_Low = results_70_NonEur$r2_low,R2_High = results_70_NonEur$r2_high)
results_70_NonEur$Train_Size <- nrow(Y_train[[1]])

results_70_EAS <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_EAS$method,R2 = results_70_EAS$r2,R2_Low = results_70_EAS$r2_low,R2_High = results_70_EAS$r2_high)
results_70_EAS$Train_Size <- nrow(Y_train[[1]])

results_70_AFR <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_AFR$method,R2 = results_70_AFR$r2,R2_Low = results_70_AFR$r2_low,R2_High = results_70_AFR$r2_high)
results_70_AFR$Train_Size <- nrow(Y_train[[1]])

results_70_SAS <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_SAS$method,R2 = results_70_SAS$r2,R2_Low = results_70_SAS$r2_low,R2_High = results_70_SAS$r2_high)
results_70_SAS$Train_Size <- nrow(Y_train[[1]])

results_70_MIX <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_MIX$method,R2 = results_70_MIX$r2,R2_Low = results_70_MIX$r2_low,R2_High = results_70_MIX$r2_high)
results_70_MIX$Train_Size <- nrow(Y_train[[1]])

results_70_UNK <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70_UNK$method,R2 = results_70_UNK$r2,R2_Low = results_70_UNK$r2_low,R2_High = results_70_UNK$r2_high)
results_70_UNK$Train_Size <- nrow(Y_train[[1]])




load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")

results_35_EUR <- NULL
results_35_NonEur <- NULL
results_35_EAS <- NULL
results_35_AFR <- NULL
results_35_SAS <- NULL
results_35_MIX <- NULL
results_35_UNK <- NULL

for(i in 1:length(Y_train)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_Eur",i,".RData"))
  sl_all_Eur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_NonEur",i,".RData"))
  sl_all_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_EAS",i,".RData"))
  sl_all_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_AFR",i,".RData"))
  sl_all_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_SAS",i,".RData"))
  sl_all_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_MIX",i,".RData"))
  sl_all_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All_UNK",i,".RData"))
  sl_all_UNK <- SL.result
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_EUR",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_EUR <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_NonEur",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_NonEur <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_EAS",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_EAS <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_AFR",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_AFR <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_SAS",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_SAS <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_MIX",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_MIX <- CT.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result_UNK",i,".RData"))
  CT.result[[1]]$method <- "CT"
  ct_UNK <- CT.result[[1]]
  
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_EUR",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_EUR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_NonEur",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_NonEur <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_EAS",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_EAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_AFR",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_AFR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_SAS",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_SAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_MIX",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_MIX <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/LASSOSUM2_result_UNK",i,".RData"))
  LASSOSUM2.result[[1]]$method <- "LASSOSUM2"
  lassosum2_UNK <- LASSOSUM2.result[[1]]
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_EUR",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_EUR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_NonEur",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_NonEur <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_EAS",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_EAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_AFR",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_AFR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_SAS",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_SAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_MIX",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_MIX <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result_UNK",i,".RData"))
  ldpred2.result[[1]]$method <- "LDPred2"
  ldpred2_UNK <- ldpred2.result[[1]]
  
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_EUR",i,".RData"))
  GeneCentric_Coding_STAARO_EUR <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_NonEUR",i,".RData"))
  GeneCentric_Coding_STAARO_NonEur <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_EAS",i,".RData"))
  GeneCentric_Coding_STAARO_EAS <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_AFR",i,".RData"))
  GeneCentric_Coding_STAARO_AFR <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_SAS",i,".RData"))
  GeneCentric_Coding_STAARO_SAS <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_MIX",i,".RData"))
  GeneCentric_Coding_STAARO_MIX <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result_UNK",i,".RData"))
  GeneCentric_Coding_STAARO_UNK <- r2.result_GeneCentric_Coding_STAARO
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_EUR",i,".RData"))
  GeneCentric_Coding_Burden_EUR <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_NonEUR",i,".RData"))
  GeneCentric_Coding_Burden_NonEur <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_EAS",i,".RData"))
  GeneCentric_Coding_Burden_EAS <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_AFR",i,".RData"))
  GeneCentric_Coding_Burden_AFR <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_SAS",i,".RData"))
  GeneCentric_Coding_Burden_SAS <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_MIX",i,".RData"))
  GeneCentric_Coding_Burden_MIX <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result_UNK",i,".RData"))
  GeneCentric_Coding_Burden_UNK <- r2.result_GeneCentric_Coding_Burden
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_EUR",i,".RData"))
  GeneCentric_Noncoding_STAARO_EUR <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_NonEUR",i,".RData"))
  GeneCentric_Noncoding_STAARO_NonEur <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_EAS",i,".RData"))
  GeneCentric_Noncoding_STAARO_EAS <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_AFR",i,".RData"))
  GeneCentric_Noncoding_STAARO_AFR <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_SAS",i,".RData"))
  GeneCentric_Noncoding_STAARO_SAS <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_MIX",i,".RData"))
  GeneCentric_Noncoding_STAARO_MIX <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result_UNK",i,".RData"))
  GeneCentric_Noncoding_STAARO_UNK <- r2.result_GeneCentric_Noncoding_STAARO
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_EUR",i,".RData"))
  GeneCentric_Noncoding_Burden_EUR <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_NonEUR",i,".RData"))
  GeneCentric_Noncoding_Burden_NonEur <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_EAS",i,".RData"))
  GeneCentric_Noncoding_Burden_EAS <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_AFR",i,".RData"))
  GeneCentric_Noncoding_Burden_AFR <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_SAS",i,".RData"))
  GeneCentric_Noncoding_Burden_SAS <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_MIX",i,".RData"))
  GeneCentric_Noncoding_Burden_MIX <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result_UNK",i,".RData"))
  GeneCentric_Noncoding_Burden_UNK <- r2.result_GeneCentric_Noncoding_Burden
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_EUR",i,".RData"))
  SlidingWindow_STAARO_EUR <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_NonEUR",i,".RData"))
  SlidingWindow_STAARO_NonEur <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_EAS",i,".RData"))
  SlidingWindow_STAARO_EAS <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_AFR",i,".RData"))
  SlidingWindow_STAARO_AFR <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_SAS",i,".RData"))
  SlidingWindow_STAARO_SAS <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_MIX",i,".RData"))
  SlidingWindow_STAARO_MIX <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result_UNK",i,".RData"))
  SlidingWindow_STAARO_UNK <- r2.result_SlidingWindow_STAARO
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_EUR",i,".RData"))
  SlidingWindow_Burden_EUR <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_NonEUR",i,".RData"))
  SlidingWindow_Burden_NonEur <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_EAS",i,".RData"))
  SlidingWindow_Burden_EAS <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_AFR",i,".RData"))
  SlidingWindow_Burden_AFR <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_SAS",i,".RData"))
  SlidingWindow_Burden_SAS <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_MIX",i,".RData"))
  SlidingWindow_Burden_MIX <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result_UNK",i,".RData"))
  SlidingWindow_Burden_UNK <- r2.result_SlidingWindow_Burden
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_EUR",i,".RData"))
  BestAll_RareVariant_STAARO_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_NonEUR",i,".RData"))
  BestAll_RareVariant_STAARO_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_EAS",i,".RData"))
  BestAll_RareVariant_STAARO_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_AFR",i,".RData"))
  BestAll_RareVariant_STAARO_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_SAS",i,".RData"))
  BestAll_RareVariant_STAARO_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_MIX",i,".RData"))
  BestAll_RareVariant_STAARO_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO_UNK",i,".RData"))
  BestAll_RareVariant_STAARO_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_EUR",i,".RData"))
  BestAll_RareVariant_Burden_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_NonEUR",i,".RData"))
  BestAll_RareVariant_Burden_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_EAS",i,".RData"))
  BestAll_RareVariant_Burden_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_AFR",i,".RData"))
  BestAll_RareVariant_Burden_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_SAS",i,".RData"))
  BestAll_RareVariant_Burden_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_MIX",i,".RData"))
  BestAll_RareVariant_Burden_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden_UNK",i,".RData"))
  BestAll_RareVariant_Burden_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_EUR",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_NonEur",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_EAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_AFR",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_SAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_MIX",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result_UNK",i,".RData"))
  SL.result$method <- "CV_plus_RV_STAARO"
  BestAll_STAARO_UNK <- SL.result
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_EUR",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_NonEur",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_EAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_AFR",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_SAS",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_MIX",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result_UNK",i,".RData"))
  SL.result$method <- "CV_plus_RV_Burden"
  BestAll_Burden_UNK <- SL.result
  
  
  results_tmp_EUR <- rbind(sl_all_Eur,ct_EUR,lassosum2_EUR,ldpred2_EUR,
                           GeneCentric_Coding_STAARO_EUR,GeneCentric_Coding_Burden_EUR,
                           GeneCentric_Noncoding_STAARO_EUR,GeneCentric_Noncoding_Burden_EUR,
                           SlidingWindow_STAARO_EUR,SlidingWindow_Burden_EUR,
                           BestAll_RareVariant_STAARO_EUR,BestAll_RareVariant_Burden_EUR,BestAll_Burden_EUR,BestAll_STAARO_EUR)
  
  results_35_EUR <- rbind(results_35_EUR,results_tmp_EUR)
  
  
  results_tmp_NonEur <- rbind(sl_all_NonEur,ct_NonEur,lassosum2_NonEur,ldpred2_NonEur,
                              GeneCentric_Coding_STAARO_NonEur,GeneCentric_Coding_Burden_NonEur,
                              GeneCentric_Noncoding_STAARO_NonEur,GeneCentric_Noncoding_Burden_NonEur,
                              SlidingWindow_STAARO_NonEur,SlidingWindow_Burden_NonEur,
                              BestAll_RareVariant_STAARO_NonEur,BestAll_RareVariant_Burden_NonEur,BestAll_Burden_NonEur,BestAll_STAARO_NonEur)
  
  results_35_NonEur <- rbind(results_35_NonEur,results_tmp_NonEur)
  
  
  results_tmp_AFR <- rbind(sl_all_AFR,ct_AFR,lassosum2_AFR,ldpred2_AFR,
                           GeneCentric_Coding_STAARO_AFR,GeneCentric_Coding_Burden_AFR,
                           GeneCentric_Noncoding_STAARO_AFR,GeneCentric_Noncoding_Burden_AFR,
                           SlidingWindow_STAARO_AFR,SlidingWindow_Burden_AFR,
                           BestAll_RareVariant_STAARO_AFR,BestAll_RareVariant_Burden_AFR,BestAll_Burden_AFR,BestAll_STAARO_AFR)
  
  results_35_AFR <- rbind(results_35_AFR,results_tmp_AFR)
  
  results_tmp_SAS <- rbind(sl_all_SAS,ct_SAS,lassosum2_SAS,ldpred2_SAS,
                           GeneCentric_Coding_STAARO_SAS,GeneCentric_Coding_Burden_SAS,
                           GeneCentric_Noncoding_STAARO_SAS,GeneCentric_Noncoding_Burden_SAS,
                           SlidingWindow_STAARO_SAS,SlidingWindow_Burden_SAS,
                           BestAll_RareVariant_STAARO_SAS,BestAll_RareVariant_Burden_SAS,BestAll_Burden_SAS,BestAll_STAARO_SAS)
  
  results_35_SAS <- rbind(results_35_SAS,results_tmp_SAS)
  
  results_tmp_EAS <- rbind(sl_all_EAS,ct_EAS,lassosum2_EAS,ldpred2_EAS,
                           GeneCentric_Coding_STAARO_EAS,GeneCentric_Coding_Burden_EAS,
                           GeneCentric_Noncoding_STAARO_EAS,GeneCentric_Noncoding_Burden_EAS,
                           SlidingWindow_STAARO_EAS,SlidingWindow_Burden_EAS,
                           BestAll_RareVariant_STAARO_EAS,BestAll_RareVariant_Burden_EAS,BestAll_Burden_EAS,BestAll_STAARO_EAS)
  
  results_35_EAS <- rbind(results_35_EAS,results_tmp_EAS)
  
  results_tmp_MIX <- rbind(sl_all_MIX,ct_MIX,lassosum2_MIX,ldpred2_MIX,
                           GeneCentric_Coding_STAARO_MIX,GeneCentric_Coding_Burden_MIX,
                           GeneCentric_Noncoding_STAARO_MIX,GeneCentric_Noncoding_Burden_MIX,
                           SlidingWindow_STAARO_MIX,SlidingWindow_Burden_MIX,
                           BestAll_RareVariant_STAARO_MIX,BestAll_RareVariant_Burden_MIX,BestAll_Burden_MIX,BestAll_STAARO_MIX)
  
  results_35_MIX <- rbind(results_35_MIX,results_tmp_MIX)
  
  results_tmp_UNK <- rbind(sl_all_UNK,ct_UNK,lassosum2_UNK,ldpred2_UNK,
                           GeneCentric_Coding_STAARO_UNK,GeneCentric_Coding_Burden_UNK,
                           GeneCentric_Noncoding_STAARO_UNK,GeneCentric_Noncoding_Burden_UNK,
                           SlidingWindow_STAARO_UNK,SlidingWindow_Burden_UNK,
                           BestAll_RareVariant_STAARO_UNK,BestAll_RareVariant_Burden_UNK,BestAll_Burden_UNK,BestAll_STAARO_UNK)
  
  results_35_UNK <- rbind(results_35_UNK,results_tmp_UNK)
  
  
  rm(list=setdiff(ls(), c("results_70_EUR","results_70_NonEur","results_70_EAS","results_70_AFR","results_70_SAS","results_70_MIX","results_70_UNK",
                          "results_35_EUR","results_35_NonEur","results_35_EAS","results_35_AFR","results_35_SAS","results_35_MIX","results_35_UNK",
                          "i","Y_train"))) 
}

scale <- rep(rep(c("Unscaled","Scaled"),each = 20*14),5)
causal_prop <- rep(c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01","Causal Prop. 0.001","Causal Prop. 0.0005"),each = 40*14)

results_35_EUR <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_EUR$method,R2 = results_35_EUR$r2,R2_Low = results_35_EUR$r2_low,R2_High = results_35_EUR$r2_high)
results_35_EUR$Train_Size <- nrow(Y_train[[1]])

results_35_NonEur <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_NonEur$method,R2 = results_35_NonEur$r2,R2_Low = results_35_NonEur$r2_low,R2_High = results_35_NonEur$r2_high)
results_35_NonEur$Train_Size <- nrow(Y_train[[1]])

results_35_EAS <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_EAS$method,R2 = results_35_EAS$r2,R2_Low = results_35_EAS$r2_low,R2_High = results_35_EAS$r2_high)
results_35_EAS$Train_Size <- nrow(Y_train[[1]])

results_35_AFR <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_AFR$method,R2 = results_35_AFR$r2,R2_Low = results_35_AFR$r2_low,R2_High = results_35_AFR$r2_high)
results_35_AFR$Train_Size <- nrow(Y_train[[1]])

results_35_SAS <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_SAS$method,R2 = results_35_SAS$r2,R2_Low = results_35_SAS$r2_low,R2_High = results_35_SAS$r2_high)
results_35_SAS$Train_Size <- nrow(Y_train[[1]])

results_35_MIX <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_MIX$method,R2 = results_35_MIX$r2,R2_Low = results_35_MIX$r2_low,R2_High = results_35_MIX$r2_high)
results_35_MIX$Train_Size <- nrow(Y_train[[1]])

results_35_UNK <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35_UNK$method,R2 = results_35_UNK$r2,R2_Low = results_35_UNK$r2_low,R2_High = results_35_UNK$r2_high)
results_35_UNK$Train_Size <- nrow(Y_train[[1]])








results_EUR <- rbind(results_35_EUR,results_70_EUR)
results_EUR$Ancestry <- "EUR"
results_NonEur <- rbind(results_35_NonEur,results_70_NonEur)
results_NonEur$Ancestry <- "NonEur"
results_EAS <- rbind(results_35_EAS,results_70_EAS)
results_EAS$Ancestry <- "EAS"
results_AFR <- rbind(results_35_AFR,results_70_AFR)
results_AFR$Ancestry <- "AFR"
results_SAS <- rbind(results_35_SAS,results_70_SAS)
results_SAS$Ancestry <- "SAS"
results_MIX <- rbind(results_35_MIX,results_70_MIX)
results_MIX$Ancestry <- "MIX"
results_UNK <- rbind(results_35_UNK,results_70_UNK)
results_UNK$Ancestry <- "UNK"


results <- rbind(results_EUR,results_NonEur,results_EAS,results_AFR,results_SAS,results_MIX,results_UNK)
results$Train_Size <- format(results$Train_Size,big.mark=",", trim=TRUE)

results$Train_Size <- paste0("n = ",results$Train_Size)

rm(list=setdiff(ls(), c("results"))) 

results <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results,mean)

results_EUR <- results[results$Ancestry == "EUR",]

rarevariant_results_EUR <- results_EUR[results_EUR$Method %in% c("GeneCentric_Coding_Burden_EUR","GeneCentric_Coding_STAARO_EUR","GeneCentric_Noncoding_Burden_EUR","GeneCentric_Noncoding_STAARO_EUR",
                                                     "SlidingWindow_Burden_EUR","SlidingWindow_STAARO_EUR","SL_All_STAARO_EUR","SL_All_Burden_EUR"),]

cv_results_EUR <- results_EUR[results_EUR$Method %in% c("CT","LASSOSUM2","LDPred2","SL_Combined_Eur"),]

overall_results_EUR <- results_EUR[results_EUR$Method %in% c("CT","LASSOSUM2","LDPred2","CV_plus_RV_STAARO","CV_plus_RV_Burden"),]

overall_results_all <- results[results$Method %in% c("CT","LASSOSUM2","LDPred2","CV_plus_RV_STAARO","CV_plus_RV_Burden"),]

overall_results_all$Method <- factor(overall_results_all$Method,levels = c("CT","LASSOSUM2","LDPred2","CV_plus_RV_Burden","CV_plus_RV_STAARO"))
overall_results_EUR$Method <- factor(overall_results_EUR$Method,levels = c("CT","LASSOSUM2","LDPred2","CV_plus_RV_Burden","CV_plus_RV_STAARO"))


####################################################### Plots

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_blank(), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

library(ggplot2)

ggplot(rarevariant_results_EUR[(rarevariant_results_EUR$Scale == "Scaled") & (rarevariant_results_EUR$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar(aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Rare Variant Methods; Scaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(rarevariant_results_EUR[(rarevariant_results_EUR$Scale == "Unscaled") & (rarevariant_results_EUR$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) +
  ggtitle("Rare Variant Methods; Unscaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()



ggplot(cv_results_EUR[(cv_results_EUR$Scale == "Scaled") & (cv_results_EUR$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Common Variant Methods; Scaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(cv_results_EUR[(cv_results_EUR$Scale == "Unscaled") & (cv_results_EUR$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Common Variant Methods; Unscaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()



ggplot(overall_results_EUR[(overall_results_EUR$Scale == "Scaled") & (overall_results_EUR$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Overall Results; Scaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(overall_results_EUR[(overall_results_EUR$Scale == "Unscaled") & (overall_results_EUR$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Overall Results; Unscaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()


ggplot(overall_results_all[(overall_results_all$Scale == "Scaled") & (overall_results_all$Train_Size == "n = 45,767") & (overall_results_all$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Overall Results; Scaled G, n = 45,767") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(overall_results_all[(overall_results_all$Scale == "Scaled") & (overall_results_all$Train_Size == "n = 91,534") & (overall_results_all$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Overall Results; Scaled G, n = 91,534") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(overall_results_all[(overall_results_all$Scale == "Unscaled") & (overall_results_all$Train_Size == "n = 45,767") & (overall_results_all$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Overall Results; Unscaled G, n = 45,767") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(overall_results_all[(overall_results_all$Scale == "Unscaled") & (overall_results_all$Train_Size == "n = 91,534") & (overall_results_all$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.01","Causal Prop. 0.0005")),]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Overall Results; Unscaled G, n = 91,534") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

