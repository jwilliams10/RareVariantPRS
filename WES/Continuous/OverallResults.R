rm(list = ls())

trait <- "BMI"

results <- NULL

for(trait in c("BMI","TC","LDL","HDL","logTG","Height")){
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EUR.RData"))
  ct_EUR <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_NonEur.RData"))
  ct_NonEur <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_AFR.RData"))
  ct_AFR <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EAS.RData"))
  ct_EAS <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_SAS.RData"))
  ct_SAS <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_MIX.RData"))
  ct_MIX <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_UNK.RData"))
  ct_UNK <- ct.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EUR.RData"))
  ldpred2_EUR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_NonEur.RData"))
  ldpred2_NonEur <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_AFR.RData"))
  ldpred2_AFR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EAS.RData"))
  ldpred2_EAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_SAS.RData"))
  ldpred2_SAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_MIX.RData"))
  ldpred2_MIX <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_UNK.RData"))
  ldpred2_UNK <- ldpred2.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EUR.RData"))
  LASSOSUM2_EUR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_NonEur.RData"))
  LASSOSUM2_NonEur <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_AFR.RData"))
  LASSOSUM2_AFR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EAS.RData"))
  LASSOSUM2_EAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_SAS.RData"))
  LASSOSUM2_SAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_MIX.RData"))
  LASSOSUM2_MIX <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_UNK.RData"))
  LASSOSUM2_UNK <- LASSOSUM2.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur.RData"))
  sl_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur.RData"))
  sl_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR.RData"))
  sl_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS.RData"))
  sl_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS.RData"))
  sl_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX.RData"))
  sl_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK.RData"))
  sl_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EUR_result.RData"))
  genecentric_coding_burden_EUR <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_NonEur_result.RData"))
  genecentric_coding_burden_NonEur <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_AFR_result.RData"))
  genecentric_coding_burden_AFR <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_EAS_result.RData"))
  genecentric_coding_burden_EAS <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_SAS_result.RData"))
  genecentric_coding_burden_SAS <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_MIX_result.RData"))
  genecentric_coding_burden_MIX <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_Burden_UNK_result.RData"))
  genecentric_coding_burden_UNK <- r2.result_GeneCentric_Coding_Burden
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EUR_result.RData"))
  genecentric_coding_STAARO_EUR <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_NonEur_result.RData"))
  genecentric_coding_STAARO_NonEur <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_AFR_result.RData"))
  genecentric_coding_STAARO_AFR <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_EAS_result.RData"))
  genecentric_coding_STAARO_EAS <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_SAS_result.RData"))
  genecentric_coding_STAARO_SAS <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_MIX_result.RData"))
  genecentric_coding_STAARO_MIX <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Coding_STAARO_UNK_result.RData"))
  genecentric_coding_STAARO_UNK <- r2.result_GeneCentric_Coding_STAARO
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EUR_result.RData"))
  GeneCentric_Noncoding_burden_EUR <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_NonEur_result.RData"))
  GeneCentric_Noncoding_burden_NonEur <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_AFR_result.RData"))
  GeneCentric_Noncoding_burden_AFR <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_EAS_result.RData"))
  GeneCentric_Noncoding_burden_EAS <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_SAS_result.RData"))
  GeneCentric_Noncoding_burden_SAS <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_MIX_result.RData"))
  GeneCentric_Noncoding_burden_MIX <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_Burden_UNK_result.RData"))
  GeneCentric_Noncoding_burden_UNK <- r2.result_GeneCentric_Noncoding_Burden
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EUR_result.RData"))
  GeneCentric_Noncoding_STAARO_EUR <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_NonEur_result.RData"))
  GeneCentric_Noncoding_STAARO_NonEur <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_AFR_result.RData"))
  GeneCentric_Noncoding_STAARO_AFR <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_EAS_result.RData"))
  GeneCentric_Noncoding_STAARO_EAS <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_SAS_result.RData"))
  GeneCentric_Noncoding_STAARO_SAS <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_MIX_result.RData"))
  GeneCentric_Noncoding_STAARO_MIX <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_GeneCentric_Noncoding_STAARO_UNK_result.RData"))
  GeneCentric_Noncoding_STAARO_UNK <- r2.result_GeneCentric_Noncoding_STAARO
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EUR_result.RData"))
  SlidingWindow_burden_EUR <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_NonEur_result.RData"))
  SlidingWindow_burden_NonEur <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_AFR_result.RData"))
  SlidingWindow_burden_AFR <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_EAS_result.RData"))
  SlidingWindow_burden_EAS <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_SAS_result.RData"))
  SlidingWindow_burden_SAS <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_MIX_result.RData"))
  SlidingWindow_burden_MIX <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_Burden_UNK_result.RData"))
  SlidingWindow_burden_UNK <- r2.result_SlidingWindow_Burden
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EUR_result.RData"))
  SlidingWindow_STAARO_EUR <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_NonEur_result.RData"))
  SlidingWindow_STAARO_NonEur <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_AFR_result.RData"))
  SlidingWindow_STAARO_AFR <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_EAS_result.RData"))
  SlidingWindow_STAARO_EAS <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_SAS_result.RData"))
  SlidingWindow_STAARO_SAS <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_MIX_result.RData"))
  SlidingWindow_STAARO_MIX <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_SlidingWindow_STAARO_UNK_result.RData"))
  SlidingWindow_STAARO_UNK <- r2.result_SlidingWindow_STAARO
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_Eur_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_NonEur_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_AFR_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_EAS_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_SAS_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_MIX_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_UNK_Burden.RData"))
  SL.result$method <- "SL_Rare_Burden"
  SL_Burden_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_Eur_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_NonEur_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_AFR_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_EAS_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_SAS_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_MIX_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_sl_result_All_UNK_STAARO.RData"))
  SL.result$method <- "SL_Rare_STAARO"
  SL_STAARO_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_EUR.RData"))
  Combined_Burden_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_NonEur.RData"))
  Combined_Burden_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_AFR.RData"))
  Combined_Burden_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_EAS.RData"))
  Combined_Burden_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_SAS.RData"))
  Combined_Burden_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_MIX.RData"))
  Combined_Burden_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_Burden_All_Result_UNK.RData"))
  Combined_Burden_UNK <- SL.result
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_EUR.RData"))
  Combined_STAARO_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_NonEur.RData"))
  Combined_STAARO_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_AFR.RData"))
  Combined_STAARO_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_EAS.RData"))
  Combined_STAARO_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_SAS.RData"))
  Combined_STAARO_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_MIX.RData"))
  Combined_STAARO_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"_STAARO_All_Result_UNK.RData"))
  Combined_STAARO_UNK <- SL.result
  
  
  
  
  results_EUR <- rbind(sl_EUR,ct_EUR,LASSOSUM2_EUR,ldpred2_EUR,genecentric_coding_burden_EUR,genecentric_coding_STAARO_EUR,GeneCentric_Noncoding_burden_EUR,
                       GeneCentric_Noncoding_STAARO_EUR,SlidingWindow_burden_EUR,SlidingWindow_STAARO_EUR,SL_Burden_EUR,SL_STAARO_EUR,Combined_Burden_EUR,Combined_STAARO_EUR)
  results_EUR$Ancestry <- "EUR"
  results_NonEur <- rbind(sl_NonEur,ct_NonEur,LASSOSUM2_NonEur,ldpred2_NonEur,genecentric_coding_burden_NonEur,genecentric_coding_STAARO_NonEur,GeneCentric_Noncoding_burden_NonEur,
                          GeneCentric_Noncoding_STAARO_NonEur,SlidingWindow_burden_NonEur,SlidingWindow_STAARO_NonEur,SL_Burden_NonEur,SL_STAARO_NonEur,Combined_Burden_NonEur,Combined_STAARO_NonEur)
  results_NonEur$Ancestry <- "NonEUR"
  results_AFR <- rbind(sl_AFR,ct_AFR,LASSOSUM2_AFR,ldpred2_AFR,genecentric_coding_burden_AFR,genecentric_coding_STAARO_AFR,GeneCentric_Noncoding_burden_AFR,
                       GeneCentric_Noncoding_STAARO_AFR,SlidingWindow_burden_AFR,SlidingWindow_STAARO_AFR,SL_Burden_AFR,SL_STAARO_AFR,Combined_Burden_AFR,Combined_STAARO_AFR)
  results_AFR$Ancestry <- "AFR"
  results_EAS <- rbind(sl_EAS,ct_EAS,LASSOSUM2_EAS,ldpred2_EAS,genecentric_coding_burden_EAS,genecentric_coding_STAARO_EAS,GeneCentric_Noncoding_burden_EAS,
                       GeneCentric_Noncoding_STAARO_EAS,SlidingWindow_burden_EAS,SlidingWindow_STAARO_EAS,SL_Burden_EAS,SL_STAARO_EAS,Combined_Burden_EAS,Combined_STAARO_EAS)
  results_EAS$Ancestry <- "EAS"
  results_SAS <- rbind(sl_SAS,ct_SAS,LASSOSUM2_SAS,ldpred2_SAS,genecentric_coding_burden_SAS,genecentric_coding_STAARO_SAS,GeneCentric_Noncoding_burden_SAS,
                       GeneCentric_Noncoding_STAARO_SAS,SlidingWindow_burden_SAS,SlidingWindow_STAARO_SAS,SL_Burden_SAS,SL_STAARO_SAS,Combined_Burden_SAS,Combined_STAARO_SAS)
  results_SAS$Ancestry <- "SAS"
  results_MIX <- rbind(sl_MIX,ct_MIX,LASSOSUM2_MIX,ldpred2_MIX,genecentric_coding_burden_MIX,genecentric_coding_STAARO_MIX,GeneCentric_Noncoding_burden_MIX,
                       GeneCentric_Noncoding_STAARO_MIX,SlidingWindow_burden_MIX,SlidingWindow_STAARO_MIX,SL_Burden_MIX,SL_STAARO_MIX,Combined_Burden_MIX,Combined_STAARO_MIX)
  results_MIX$Ancestry <- "MIX"
  results_UNK <- rbind(sl_UNK,ct_UNK,LASSOSUM2_UNK,ldpred2_UNK,genecentric_coding_burden_UNK,genecentric_coding_STAARO_UNK,GeneCentric_Noncoding_burden_UNK,
                       GeneCentric_Noncoding_STAARO_UNK,SlidingWindow_burden_UNK,SlidingWindow_STAARO_UNK,SL_Burden_UNK,SL_STAARO_UNK,Combined_Burden_UNK,Combined_STAARO_UNK)
  results_UNK$Ancestry <- "UNK"
  
  results_tmp <- rbind(results_EUR,results_NonEur,results_AFR,results_EAS,results_SAS,results_MIX,results_UNK)
  results_tmp$Trait <- trait
   
  results <- rbind(results,results_tmp)
}

library(stringr)

results$method[str_detect(results$method,"CT")] <- "CT"
results$method[str_detect(results$method,"SL_Combined")] <- "SL_Combined"
results$method[str_detect(results$method,"LDPred2")] <- "LDPred2"
results$method[str_detect(results$method,"LASSOSUM2")] <- "LASSOSUM2"
results$method[str_detect(results$method,"CV_plus_RV_STAARO")] <- "CV_plus_RV_STAARO"
results$method[str_detect(results$method,"CV_plus_RV_Burden")] <- "CV_plus_RV_Burden"
results$method[str_detect(results$method,"SlidingWindow_STAARO")] <- "SlidingWindow_STAARO"
results$method[str_detect(results$method,"SlidingWindow_Burden")] <- "SlidingWindow_Burden"
results$method[str_detect(results$method,"GeneCentric_Coding_STAARO")] <- "GeneCentric_Coding_STAARO"
results$method[str_detect(results$method,"GeneCentric_Coding_Burden")] <- "GeneCentric_Coding_Burden"
results$method[str_detect(results$method,"GeneCentric_Noncoding_STAARO")] <- "GeneCentric_Noncoding_STAARO"
results$method[str_detect(results$method,"GeneCentric_Noncoding_Burden")] <- "GeneCentric_Noncoding_Burden"

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

common <- results[results$method %in% c("CT","LDPred2","LASSOSUM2","SL_Combined"),]
common$method <- factor(common$method,levels = c("CT","LDPred2","LASSOSUM2","SL_Combined"))

ggplot(common) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Common Variants") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

rare <- results[results$method %in% c("SlidingWindow_STAARO","SlidingWindow_Burden","GeneCentric_Coding_STAARO","GeneCentric_Coding_Burden","GeneCentric_Noncoding_STAARO","GeneCentric_Noncoding_Burden","SL_Rare_Burden","SL_Rare_STAARO"),]
rare$method <- factor(rare$method,levels = c("GeneCentric_Coding_STAARO","GeneCentric_Coding_Burden","GeneCentric_Noncoding_STAARO","GeneCentric_Noncoding_Burden","SlidingWindow_Burden","SlidingWindow_STAARO","SL_Rare_Burden","SL_Rare_STAARO"))

ggplot(rare) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Rare Variants") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

paper <- results[results$method %in% c("CT","LDPred2","LASSOSUM2","CV_plus_RV_STAARO"),]
paper$method <- factor(paper$method,levels = c("CT","LDPred2","LASSOSUM2","CV_plus_RV_STAARO"))

ggplot(paper[!(paper$Ancestry %in% c("UNK","NonEUR")),]) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Paper Plot") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

# Max 5, no error bars.

ggplot(paper[paper$Ancestry == "EUR" & paper$Trait == "HDL",]) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
#  geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
#  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Paper Plot") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()
