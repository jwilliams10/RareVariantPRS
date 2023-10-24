rm(list = ls())

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_RareVariants_PRS/final_coefs_All_STAARO.RData")

Train_Effect_Sizes_All_GeneCentric_Coding <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Train_Effect_Sizes_All.csv")

i <- 1

betas_GeneCentricCoding <- NULL

for(i in 1:length(thresholds)){
  beta_i <- vector()
  beta_i[Train_Effect_Sizes_All_GeneCentric_Coding$STAAR_O <= thresholds[i]] <- Train_Effect_Sizes_All_GeneCentric_Coding$Burden_Est[Train_Effect_Sizes_All_GeneCentric_Coding$STAAR_O <= thresholds[i]]
  beta_i[is.na(beta_i)] <- 0
  betas_GeneCentricCoding <- cbind(betas_GeneCentricCoding,beta_i)
}
colnames(betas_GeneCentricCoding) <- c(paste0("GeneCentric_Coding_PRS_Threshold_",1:15))

Train_Effect_Sizes_All_GeneCentric_Noncoding <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Train_Effect_Sizes_All.csv")

i <- 1

betas_GeneCentricNoncoding <- NULL

for(i in 1:length(thresholds)){
  beta_i <- vector()
  beta_i[Train_Effect_Sizes_All_GeneCentric_Noncoding$STAAR_O <= thresholds[i]] <- Train_Effect_Sizes_All_GeneCentric_Noncoding$Burden_Est[Train_Effect_Sizes_All_GeneCentric_Noncoding$STAAR_O <= thresholds[i]]
  beta_i[is.na(beta_i)] <- 0
  betas_GeneCentricNoncoding <- cbind(betas_GeneCentricNoncoding,beta_i)
}
colnames(betas_GeneCentricNoncoding) <- c(paste0("GeneCentric_Noncoding_PRS_Threshold_",1:15))

Train_Effect_Sizes_All_SlidingWindow <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Train_Effect_Sizes_All.csv")

i <- 1

betas_SlidingWindow <- NULL

for(i in 1:length(thresholds)){
  beta_i <- vector()
  beta_i[Train_Effect_Sizes_All_SlidingWindow$STAAR_O <= thresholds[i]] <- Train_Effect_Sizes_All_SlidingWindow$Burden_Est[Train_Effect_Sizes_All_SlidingWindow$STAAR_O <= thresholds[i]]
  beta_i[is.na(beta_i)] <- 0
  betas_SlidingWindow <- cbind(betas_SlidingWindow,beta_i)
}
colnames(betas_SlidingWindow) <- c(paste0("SlidingWindow_PRS_Threshold_",1:15))

best_beta_genecentric_coding <- betas_GeneCentricCoding[,colnames(betas_GeneCentricCoding) %in% names(final_coefs)] %*% matrix(final_coefs[names(final_coefs)%in% colnames(betas_GeneCentricCoding)],ncol = 1)
best_beta_genecentric_noncoding <- betas_GeneCentricNoncoding[,colnames(betas_GeneCentricNoncoding) %in% names(final_coefs)] %*% matrix(final_coefs[names(final_coefs)%in% colnames(betas_GeneCentricNoncoding)],ncol = 1)
best_beta_slidingwindow <- betas_SlidingWindow[,colnames(betas_SlidingWindow) %in% names(final_coefs)] %*% matrix(final_coefs[names(final_coefs)%in% colnames(betas_SlidingWindow)],ncol = 1)

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Coefficients_STAARO.RData")

best_beta_genecentric_coding <- best_beta_genecentric_coding*effects[2]
best_beta_genecentric_noncoding <- best_beta_genecentric_noncoding*effects[2]
best_beta_slidingwindow <- best_beta_slidingwindow*effects[2]

Train_Effect_Sizes_All_GeneCentric_Coding$Final_Beta <- best_beta_genecentric_coding
Train_Effect_Sizes_All_GeneCentric_Noncoding$Final_Beta <- best_beta_genecentric_noncoding
Train_Effect_Sizes_All_SlidingWindow$Final_Beta <- best_beta_slidingwindow

Train_Effect_Sizes_All_GeneCentric_Coding <- Train_Effect_Sizes_All_GeneCentric_Coding[,c("Gene","Chr","Category","Final_Beta")]
Train_Effect_Sizes_All_GeneCentric_Noncoding <- Train_Effect_Sizes_All_GeneCentric_Noncoding[,c("Gene","Chr","Category","Final_Beta")]
Train_Effect_Sizes_All_SlidingWindow <- Train_Effect_Sizes_All_SlidingWindow[,c("Chr","Start","End","Final_Beta")]

write.csv(Train_Effect_Sizes_All_SlidingWindow,row.names = FALSE,"/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Final_Coefficients.csv")
write.csv(Train_Effect_Sizes_All_GeneCentric_Coding,row.names = FALSE,"/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Final_Coefficients.csv")
write.csv(Train_Effect_Sizes_All_GeneCentric_Noncoding,row.names = FALSE,"/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Final_Coefficients.csv")
