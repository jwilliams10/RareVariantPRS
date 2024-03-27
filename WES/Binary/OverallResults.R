rm(list = ls())

trait <- "Asthma"

results <- NULL

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_EUR.RData"))
  ct_EUR <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_NonEur.RData"))
  ct_NonEur <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_AFR.RData"))
  ct_AFR <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_EAS.RData"))
  ct_EAS <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_SAS.RData"))
  ct_SAS <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_MIX.RData"))
  ct_MIX <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_CT_result_UNK.RData"))
  ct_UNK <- ct.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_EUR.RData"))
  ldpred2_EUR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_NonEur.RData"))
  ldpred2_NonEur <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_AFR.RData"))
  ldpred2_AFR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_EAS.RData"))
  ldpred2_EAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_SAS.RData"))
  ldpred2_SAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_MIX.RData"))
  ldpred2_MIX <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_UNK.RData"))
  ldpred2_UNK <- ldpred2.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EUR.RData"))
  LASSOSUM2_EUR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_NonEur.RData"))
  LASSOSUM2_NonEur <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_AFR.RData"))
  LASSOSUM2_AFR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EAS.RData"))
  LASSOSUM2_EAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_SAS.RData"))
  LASSOSUM2_SAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_MIX.RData"))
  LASSOSUM2_MIX <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_UNK.RData"))
  LASSOSUM2_UNK <- LASSOSUM2.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur.RData"))
  sl_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur.RData"))
  sl_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR.RData"))
  sl_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS.RData"))
  sl_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS.RData"))
  sl_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX.RData"))
  sl_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK.RData"))
  sl_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EUR_result.RData"))
  genecentric_coding_burden_EUR <- AUC.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_NonEur_result.RData"))
  genecentric_coding_burden_NonEur <- AUC.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_AFR_result.RData"))
  genecentric_coding_burden_AFR <- AUC.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EAS_result.RData"))
  genecentric_coding_burden_EAS <- AUC.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_SAS_result.RData"))
  genecentric_coding_burden_SAS <- AUC.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_MIX_result.RData"))
  genecentric_coding_burden_MIX <- AUC.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_UNK_result.RData"))
  genecentric_coding_burden_UNK <- AUC.result_GeneCentric_Coding_Burden
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EUR_result.RData"))
  genecentric_coding_STAARO_EUR <- AUC.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_NonEur_result.RData"))
  genecentric_coding_STAARO_NonEur <- AUC.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_AFR_result.RData"))
  genecentric_coding_STAARO_AFR <- AUC.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EAS_result.RData"))
  genecentric_coding_STAARO_EAS <- AUC.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_SAS_result.RData"))
  genecentric_coding_STAARO_SAS <- AUC.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_MIX_result.RData"))
  genecentric_coding_STAARO_MIX <- AUC.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_UNK_result.RData"))
  genecentric_coding_STAARO_UNK <- AUC.result_GeneCentric_Coding_STAARO
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EUR_result.RData"))
  GeneCentric_Noncoding_burden_EUR <- AUC.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_NonEur_result.RData"))
  GeneCentric_Noncoding_burden_NonEur <- AUC.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_AFR_result.RData"))
  GeneCentric_Noncoding_burden_AFR <- AUC.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EAS_result.RData"))
  GeneCentric_Noncoding_burden_EAS <- AUC.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_SAS_result.RData"))
  GeneCentric_Noncoding_burden_SAS <- AUC.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_MIX_result.RData"))
  GeneCentric_Noncoding_burden_MIX <- AUC.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_UNK_result.RData"))
  GeneCentric_Noncoding_burden_UNK <- AUC.result_GeneCentric_Noncoding_Burden
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EUR_result.RData"))
  GeneCentric_Noncoding_STAARO_EUR <- AUC.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_NonEur_result.RData"))
  GeneCentric_Noncoding_STAARO_NonEur <- AUC.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_AFR_result.RData"))
  GeneCentric_Noncoding_STAARO_AFR <- AUC.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EAS_result.RData"))
  GeneCentric_Noncoding_STAARO_EAS <- AUC.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_SAS_result.RData"))
  GeneCentric_Noncoding_STAARO_SAS <- AUC.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_MIX_result.RData"))
  GeneCentric_Noncoding_STAARO_MIX <- AUC.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_UNK_result.RData"))
  GeneCentric_Noncoding_STAARO_UNK <- AUC.result_GeneCentric_Noncoding_STAARO
  
  
  
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EUR_result.RData"))
  # SlidingWindow_burden_EUR <- AUC.result_SlidingWindow_Burden
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_NonEur_result.RData"))
  # SlidingWindow_burden_NonEur <- AUC.result_SlidingWindow_Burden
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_AFR_result.RData"))
  # SlidingWindow_burden_AFR <- AUC.result_SlidingWindow_Burden
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EAS_result.RData"))
  # SlidingWindow_burden_EAS <- AUC.result_SlidingWindow_Burden
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_SAS_result.RData"))
  # SlidingWindow_burden_SAS <- AUC.result_SlidingWindow_Burden
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_MIX_result.RData"))
  # SlidingWindow_burden_MIX <- AUC.result_SlidingWindow_Burden
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_UNK_result.RData"))
  # SlidingWindow_burden_UNK <- AUC.result_SlidingWindow_Burden
  
  
  
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EUR_result.RData"))
  # SlidingWindow_STAARO_EUR <- AUC.result_SlidingWindow_STAARO
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_NonEur_result.RData"))
  # SlidingWindow_STAARO_NonEur <- AUC.result_SlidingWindow_STAARO
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_AFR_result.RData"))
  # SlidingWindow_STAARO_AFR <- AUC.result_SlidingWindow_STAARO
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EAS_result.RData"))
  # SlidingWindow_STAARO_EAS <- AUC.result_SlidingWindow_STAARO
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_SAS_result.RData"))
  # SlidingWindow_STAARO_SAS <- AUC.result_SlidingWindow_STAARO
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_MIX_result.RData"))
  # SlidingWindow_STAARO_MIX <- AUC.result_SlidingWindow_STAARO
  # 
  # load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_UNK_result.RData"))
  # SlidingWindow_STAARO_UNK <- AUC.result_SlidingWindow_STAARO
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_Eur_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_EUR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_NonEur_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_NonEur <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_AFR_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_AFR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_EAS_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_EAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_SAS_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_SAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_MIX_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_MIX <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_UNK_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_UNK <- SL.result



  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_Eur_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_EUR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_NonEur_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_NonEur <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_AFR_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_AFR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_EAS_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_EAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_SAS_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_SAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_MIX_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_MIX <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_UNK_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_UNK <- SL.result



  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_EUR.RData"))
  Combined_Burden_EUR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_NonEur.RData"))
  Combined_Burden_NonEur <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_AFR.RData"))
  Combined_Burden_AFR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_EAS.RData"))
  Combined_Burden_EAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_SAS.RData"))
  Combined_Burden_SAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_MIX.RData"))
  Combined_Burden_MIX <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_UNK.RData"))
  Combined_Burden_UNK <- SL.result



  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_EUR.RData"))
  Combined_STAARO_EUR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_NonEur.RData"))
  Combined_STAARO_NonEur <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_AFR.RData"))
  Combined_STAARO_AFR <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_EAS.RData"))
  Combined_STAARO_EAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_SAS.RData"))
  Combined_STAARO_SAS <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_MIX.RData"))
  Combined_STAARO_MIX <- SL.result

  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_UNK.RData"))
  Combined_STAARO_UNK <- SL.result
  
  
  
  
  results_EUR <- rbind(sl_EUR,ct_EUR,LASSOSUM2_EUR,ldpred2_EUR,genecentric_coding_burden_EUR,genecentric_coding_STAARO_EUR,GeneCentric_Noncoding_burden_EUR,
                       GeneCentric_Noncoding_STAARO_EUR,SL_Burden_EUR,SL_STAARO_EUR,Combined_Burden_EUR,Combined_STAARO_EUR)
  results_EUR$Ancestry <- "EUR"
  results_NonEur <- rbind(sl_NonEur,ct_NonEur,LASSOSUM2_NonEur,ldpred2_NonEur,genecentric_coding_burden_NonEur,genecentric_coding_STAARO_NonEur,GeneCentric_Noncoding_burden_NonEur,
                          GeneCentric_Noncoding_STAARO_NonEur,SL_Burden_NonEur,SL_STAARO_NonEur,Combined_Burden_NonEur,Combined_STAARO_NonEur)
  results_NonEur$Ancestry <- "NonEUR"
  results_AFR <- rbind(sl_AFR,ct_AFR,LASSOSUM2_AFR,ldpred2_AFR,genecentric_coding_burden_AFR,genecentric_coding_STAARO_AFR,GeneCentric_Noncoding_burden_AFR,
                       GeneCentric_Noncoding_STAARO_AFR,SL_Burden_AFR,SL_STAARO_AFR,Combined_Burden_AFR,Combined_STAARO_AFR)
  results_AFR$Ancestry <- "AFR"
  results_EAS <- rbind(sl_EAS,ct_EAS,LASSOSUM2_EAS,ldpred2_EAS,genecentric_coding_burden_EAS,genecentric_coding_STAARO_EAS,GeneCentric_Noncoding_burden_EAS,
                       GeneCentric_Noncoding_STAARO_EAS,SL_Burden_EAS,SL_STAARO_EAS,Combined_Burden_EAS,Combined_STAARO_EAS)
  results_EAS$Ancestry <- "EAS"
  results_SAS <- rbind(sl_SAS,ct_SAS,LASSOSUM2_SAS,ldpred2_SAS,genecentric_coding_burden_SAS,genecentric_coding_STAARO_SAS,GeneCentric_Noncoding_burden_SAS,
                       GeneCentric_Noncoding_STAARO_SAS,SL_Burden_SAS,SL_STAARO_SAS,Combined_Burden_SAS,Combined_STAARO_SAS)
  results_SAS$Ancestry <- "SAS"
  results_MIX <- rbind(sl_MIX,ct_MIX,LASSOSUM2_MIX,ldpred2_MIX,genecentric_coding_burden_MIX,genecentric_coding_STAARO_MIX,GeneCentric_Noncoding_burden_MIX,
                       GeneCentric_Noncoding_STAARO_MIX,SL_Burden_MIX,SL_STAARO_MIX,Combined_Burden_MIX,Combined_STAARO_MIX)
  results_MIX$Ancestry <- "MIX"
  results_UNK <- rbind(sl_UNK,ct_UNK,LASSOSUM2_UNK,ldpred2_UNK,genecentric_coding_burden_UNK,genecentric_coding_STAARO_UNK,GeneCentric_Noncoding_burden_UNK,
                       GeneCentric_Noncoding_STAARO_UNK,SL_Burden_UNK,SL_STAARO_UNK,Combined_Burden_UNK,Combined_STAARO_UNK)
  results_UNK$Ancestry <- "UNK"
  
  results_tmp <- rbind(results_EUR,results_NonEur,results_AFR,results_EAS,results_SAS,results_MIX,results_UNK)
  results_tmp$Trait <- trait
  
  results <- rbind(results,results_tmp)
}

library(stringr)

results$Method <- results$method
results <- subset(results,select = -method)

results$Method[str_detect(results$Method,"CT")] <- "CT"
results$Method[str_detect(results$Method,"SL_Combined")] <- "SL_Combined"
results$Method[str_detect(results$Method,"LDPred2")] <- "LDPred2"
results$Method[str_detect(results$Method,"LASSOSUM2")] <- "LASSOSUM2"
results$Method[str_detect(results$Method,"CV_plus_RV_STAARO")] <- "CV_plus_RV_STAARB"
results$Method[str_detect(results$Method,"CV_plus_RV_Burden")] <- "CV_plus_RV_Burden"
results$Method[str_detect(results$Method,"SlidingWindow_STAARO")] <- "SlidingWindow_STAARB"
results$Method[str_detect(results$Method,"SlidingWindow_Burden")] <- "SlidingWindow_Burden"
results$Method[str_detect(results$Method,"GeneCentric_Coding_STAARO")] <- "GeneCentric_Coding_STAARB"
results$Method[str_detect(results$Method,"GeneCentric_Coding_Burden")] <- "GeneCentric_Coding_Burden"
results$Method[str_detect(results$Method,"GeneCentric_Noncoding_STAARO")] <- "GeneCentric_Noncoding_STAARB"
results$Method[str_detect(results$Method,"GeneCentric_Noncoding_Burden")] <- "GeneCentric_Noncoding_Burden"

results$Method[str_detect(results$Method,"SL_Rare_STAARO")] <- "SL_Rare_STAARB"

rm(list=setdiff(ls(), "results"))

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

common <- results[results$Method %in% c("CT","LDPred2","LASSOSUM2","SL_Combined"),]
common$Method <- factor(common$Method,levels = c("CT","LDPred2","LASSOSUM2","SL_Combined"))

ggplot(common[!(common$Ancestry %in% c("EAS","UNK","MIX")),]) +
  geom_bar(aes(x=Method, y=AUC,fill=Method), stat="identity", alpha=0.7) +
#  geom_errorbar( aes(x=Method, ymin=AUC_low, ymax=AUC_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("WES Common Variants") + 
  coord_cartesian(ylim=c(0.45,0.8)) + 
  ylab("AUC") + 
  theme_Publication() + 
  scale_fill_Publication()

rare <- results[results$Method %in% c("SlidingWindow_STAARO","SlidingWindow_Burden","GeneCentric_Coding_STAARB","GeneCentric_Coding_Burden","GeneCentric_Noncoding_STAARB","GeneCentric_Noncoding_Burden","SL_Rare_Burden","SL_Rare_STAARB"),]
rare$Method <- factor(rare$Method,levels = c("GeneCentric_Coding_STAARB","GeneCentric_Coding_Burden","GeneCentric_Noncoding_STAARB","GeneCentric_Noncoding_Burden","SL_Rare_Burden","SL_Rare_STAARB"))

ggplot(rare[!(rare$Ancestry %in% c("EAS","UNK","MIX")),]) +
  geom_bar(aes(x=Method, y=AUC,fill=Method), stat="identity", alpha=0.7) +
#  geom_errorbar( aes(x=Method, ymin=AUC_low, ymax=AUC_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("WES Rare Variants") + 
  ylab("AUC") + 
  coord_cartesian(ylim=c(0.45,0.7)) + 
  theme_Publication() + 
  scale_fill_Publication()

paper <- results[results$Method %in% c("CT","LDPred2","LASSOSUM2","CV_plus_RV_STAARB","CV_plus_RV_Burden"),]
paper$Method <- factor(paper$Method,levels = c("CT","LDPred2","LASSOSUM2","CV_plus_RV_Burden","CV_plus_RV_STAARB"))

ggplot(paper[!(paper$Ancestry %in% c("EAS","UNK","MIX")),]) +
  geom_bar(aes(x=Method, y=AUC,fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=AUC_low, ymax=AUC_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("WES Overall") + 
  ylab("AUC") + 
  coord_cartesian(ylim=c(0.45,0.8)) + 
  theme_Publication() + 
  scale_fill_Publication()
