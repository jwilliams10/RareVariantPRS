rm(list = ls())

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/sl_result_All.RData")
sl_all <- SL.result
sl_all$method <- "SL_All"

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/sl_result_BestThree.RData")
sl_bestthree <- SL.result
sl_bestthree$method <- "SL_BestThree"

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/CT/CT_result.RData")
ct <- ct.result[[1]]

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2_result.RData")
lassosum2 <- lassosum2.result[[1]]

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_result.RData")
ldpred2 <- ldpred2.result[[1]]

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/GeneCentric_Coding_STAARO_result.RData")
GeneCentric_Coding_STAARO <- r2.result_GeneCentric_Coding_STAARO

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/GeneCentric_Coding_Burden_result.RData")
GeneCentric_Coding_Burden <- r2.result_GeneCentric_Coding_Burden

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/GeneCentric_Noncoding_STAARO_result.RData")
GeneCentric_Noncoding_STAARO <- r2.result_GeneCentric_Noncoding_STAARO

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/GeneCentric_Noncoding_Burden_result.RData")
GeneCentric_Noncoding_Burden <- r2.result_GeneCentric_Noncoding_Burden

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/SlidingWindow_STAARO_result.RData")
SlidingWindow_STAARO <- r2.result_SlidingWindow_STAARO

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/SlidingWindow_Burden_result.RData")
SlidingWindow_Burden <- r2.result_SlidingWindow_Burden


load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/sl_result_All_STAARO.RData")
BestAll_RareVariant_STAARO <- SL.result
BestAll_RareVariant_STAARO$method <- "SL_All_STAARO"

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/sl_result_All_Burden.RData")
BestAll_RareVariant_Burden <- SL.result
BestAll_RareVariant_Burden$method <- "SL_All_Burden"

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/sl_result_BestThree_STAARO.RData")
BestThree_RareVariant_STAARO <- SL.result
BestThree_RareVariant_STAARO$method <- "SL_BestThree_STAARO"

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/sl_result_BestThree_Burden.RData")
BestThree_RareVariant_Burden <- SL.result
BestThree_RareVariant_Burden$method <- "SL_BestThree_Burden"



BestThree_Burden <- data.frame(method = "CV_RV_Burden_BestThree", r2 = 0.1702937,r2_low = NA, r2_high = NA)
BestThree_STAARO <- data.frame(method = "CV_RV_STAARO_BestThree", r2 = 0.1705838,r2_low = NA, r2_high = NA)

BestAll_Burden <- data.frame(method = "CV_RV_Burden_All", r2 = 0.1704095,r2_low = NA, r2_high = NA)
BestAll_STAARO <- data.frame(method = "CV_RV_STAARO_All", r2 = 0.1702076,r2_low = NA, r2_high = NA)


results <- rbind(sl_all,sl_bestthree,ct,lassosum2,ldpred2,
                 GeneCentric_Coding_STAARO,GeneCentric_Coding_Burden,
                 GeneCentric_Noncoding_STAARO,GeneCentric_Noncoding_Burden,
                 SlidingWindow_STAARO,SlidingWindow_Burden,
                 BestAll_RareVariant_STAARO,BestAll_RareVariant_Burden,
                 BestThree_RareVariant_STAARO,BestThree_RareVariant_Burden,
                 BestThree_Burden,BestThree_STAARO,BestAll_Burden,BestAll_STAARO)

# results[,2:4] <- results[,2:4]*100

rm(list=setdiff(ls(), "results"))
