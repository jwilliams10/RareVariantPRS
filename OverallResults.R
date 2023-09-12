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

results <- rbind(sl_all,sl_bestthree,ct,lassosum2,ldpred2)

rm(list=setdiff(ls(), "results"))
