rm(list = ls())

library(dplyr)
library("xlsx")

final_results_WES <- NULL

full_results <- NULL
full_results_Boot <- NULL

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  CT_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_Bootstraps.csv"))
  CT_Boot_Results$Method <- "CT"
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  LDPred2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_Bootstraps.csv"))
  LDPred2_Boot_Results$Method <- "LDpred2"
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  LASSOSUM2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_Bootstraps.csv"))
  LASSOSUM2_Boot_Results$Method <- "Lassosum2"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_CV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/CV_",trait,"_Bootstraps.csv"))
  RICE_CV_Boot_Results$Method <- "RICE-CV"
  colnames(RICE_CV_Boot_Results) <- colnames(CT_Boot_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  RICE_RV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/RV_",trait,"_Bootstraps.csv"))
  RICE_RV_Boot_Results$Method <- "RICE-RV"
  colnames(RICE_RV_Boot_Results) <- colnames(CT_Boot_Results)
  full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,RICE_CV_Results,RICE_RV_Results))
  full_results_Boot <- rbind(full_results_Boot,rbind(CT_Boot_Results,LDPred2_Boot_Results,LASSOSUM2_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_Boot$trait <- factor(full_results_Boot$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_SAS_boot_Lower,CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_SAS_boot_Upper,CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper),
                    R2_raw_Lower_95 = c(CI_95$R2_raw_EUR_boot_Lower,CI_95$R2_raw_SAS_boot_Lower,CI_95$R2_raw_AFR_boot_Lower,CI_95$R2_raw_AMR_boot_Lower),
                    R2_raw_Upper_95 = c(CI_95$R2_raw_EUR_boot_Upper,CI_95$R2_raw_SAS_boot_Upper,CI_95$R2_raw_AFR_boot_Upper,CI_95$R2_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower,CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper,CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper),
                    R2_adjusted_Lower_95 = c(CI_95$R2_adjusted_EUR_boot_Lower,CI_95$R2_adjusted_SAS_boot_Lower,CI_95$R2_adjusted_AFR_boot_Lower,CI_95$R2_adjusted_AMR_boot_Lower),
                    R2_adjusted_Upper_95 = c(CI_95$R2_adjusted_EUR_boot_Upper,CI_95$R2_adjusted_SAS_boot_Upper,CI_95$R2_adjusted_AFR_boot_Upper,CI_95$R2_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_adjusted","R2_adjusted_Lower_95","R2_adjusted_Upper_95")]
colnames(full_results) <- c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_AUC_adjusted","R2_AUC_adjusted_Lower_95","R2_AUC_adjusted_Upper_95")

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_WES <- rbind(final_results_WES,full_results)

full_results <- NULL
full_results_Boot <- NULL

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  CT_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_Bootstraps.csv"))
  CT_Boot_Results$Method <- "CT"
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  LDPred2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_Bootstraps.csv"))
  LDPred2_Boot_Results$Method <- "LDpred2"
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  LASSOSUM2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_Bootstraps.csv"))
  LASSOSUM2_Boot_Results$Method <- "Lassosum2"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_CV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/CV_",trait,"_Bootstraps.csv"))
  RICE_CV_Boot_Results$Method <- "RICE-CV"
  colnames(RICE_CV_Boot_Results) <- colnames(CT_Boot_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  RICE_RV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Common_plus_RareVariants/RV_",trait,"_Bootstraps.csv"))
  RICE_RV_Boot_Results$Method <- "RICE-RV"
  colnames(RICE_RV_Boot_Results) <- colnames(CT_Boot_Results)
  full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,RICE_CV_Results,RICE_RV_Results))
  full_results_Boot <- rbind(full_results_Boot,rbind(CT_Boot_Results,LDPred2_Boot_Results,LASSOSUM2_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait <- factor(full_results$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_SAS_boot_Lower,CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_SAS_boot_Upper,CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper),
                    AUC_raw_Lower_95 = c(CI_95$AUC_raw_EUR_boot_Lower,CI_95$AUC_raw_SAS_boot_Lower,CI_95$AUC_raw_AFR_boot_Lower,CI_95$AUC_raw_AMR_boot_Lower),
                    AUC_raw_Upper_95 = c(CI_95$AUC_raw_EUR_boot_Upper,CI_95$AUC_raw_SAS_boot_Upper,CI_95$AUC_raw_AFR_boot_Upper,CI_95$AUC_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower,CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper,CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper),
                    AUC_adjusted_Lower_95 = c(CI_95$AUC_adjusted_EUR_boot_Lower,CI_95$AUC_adjusted_SAS_boot_Lower,CI_95$AUC_adjusted_AFR_boot_Lower,CI_95$AUC_adjusted_AMR_boot_Lower),
                    AUC_adjusted_Upper_95 = c(CI_95$AUC_adjusted_EUR_boot_Upper,CI_95$AUC_adjusted_SAS_boot_Upper,CI_95$AUC_adjusted_AFR_boot_Upper,CI_95$AUC_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","AUC_adjusted","AUC_adjusted_Lower_95","AUC_adjusted_Upper_95")]
colnames(full_results) <- c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_AUC_adjusted","R2_AUC_adjusted_Lower_95","R2_AUC_adjusted_Upper_95")

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_WES <- rbind(final_results_WES,full_results)

write.xlsx(final_results_WES, file = "Metrics_Table.xlsx",sheetName = "UKB WES", append = TRUE,row.names = FALSE)


rm(list = ls())

final_results_Imputed <- NULL

full_results <- NULL
full_results_Boot <- NULL

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  CT_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Bootstraps.csv"))
  CT_Boot_Results$Method <- "CT"
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  LDPred2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Bootstraps.csv"))
  LDPred2_Boot_Results$Method <- "LDpred2"
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  LASSOSUM2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Bootstraps.csv"))
  LASSOSUM2_Boot_Results$Method <- "Lassosum2"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_CV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"_Bootstraps.csv"))
  RICE_CV_Boot_Results$Method <- "RICE-CV"
  colnames(RICE_CV_Boot_Results) <- colnames(CT_Boot_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  RICE_RV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"_Bootstraps.csv"))
  RICE_RV_Boot_Results$Method <- "RICE-RV"
  colnames(RICE_RV_Boot_Results) <- colnames(CT_Boot_Results)
  full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,RICE_CV_Results,RICE_RV_Results))
  full_results_Boot <- rbind(full_results_Boot,rbind(CT_Boot_Results,LDPred2_Boot_Results,LASSOSUM2_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_Boot$trait <- factor(full_results_Boot$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_SAS_boot_Lower,CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_SAS_boot_Upper,CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper),
                    R2_raw_Lower_95 = c(CI_95$R2_raw_EUR_boot_Lower,CI_95$R2_raw_SAS_boot_Lower,CI_95$R2_raw_AFR_boot_Lower,CI_95$R2_raw_AMR_boot_Lower),
                    R2_raw_Upper_95 = c(CI_95$R2_raw_EUR_boot_Upper,CI_95$R2_raw_SAS_boot_Upper,CI_95$R2_raw_AFR_boot_Upper,CI_95$R2_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower,CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper,CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper),
                    R2_adjusted_Lower_95 = c(CI_95$R2_adjusted_EUR_boot_Lower,CI_95$R2_adjusted_SAS_boot_Lower,CI_95$R2_adjusted_AFR_boot_Lower,CI_95$R2_adjusted_AMR_boot_Lower),
                    R2_adjusted_Upper_95 = c(CI_95$R2_adjusted_EUR_boot_Upper,CI_95$R2_adjusted_SAS_boot_Upper,CI_95$R2_adjusted_AFR_boot_Upper,CI_95$R2_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_adjusted","R2_adjusted_Lower_95","R2_adjusted_Upper_95")]
colnames(full_results) <- c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_AUC_adjusted","R2_AUC_adjusted_Lower_95","R2_AUC_adjusted_Upper_95")

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_Imputed <- rbind(final_results_Imputed,full_results)

full_results <- NULL
full_results_Boot <- NULL

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  CT_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Bootstraps.csv"))
  CT_Boot_Results$Method <- "CT"
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  LDPred2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Bootstraps.csv"))
  LDPred2_Boot_Results$Method <- "LDpred2"
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  LASSOSUM2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Bootstraps.csv"))
  LASSOSUM2_Boot_Results$Method <- "Lassosum2"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_CV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"_Bootstraps.csv"))
  RICE_CV_Boot_Results$Method <- "RICE-CV"
  colnames(RICE_CV_Boot_Results) <- colnames(CT_Boot_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  RICE_RV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"_Bootstraps.csv"))
  RICE_RV_Boot_Results$Method <- "RICE-RV"
  colnames(RICE_RV_Boot_Results) <- colnames(CT_Boot_Results)
  full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,RICE_CV_Results,RICE_RV_Results))
  full_results_Boot <- rbind(full_results_Boot,rbind(CT_Boot_Results,LDPred2_Boot_Results,LASSOSUM2_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait <- factor(full_results$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_SAS_boot_Lower,CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_SAS_boot_Upper,CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper),
                    AUC_raw_Lower_95 = c(CI_95$AUC_raw_EUR_boot_Lower,CI_95$AUC_raw_SAS_boot_Lower,CI_95$AUC_raw_AFR_boot_Lower,CI_95$AUC_raw_AMR_boot_Lower),
                    AUC_raw_Upper_95 = c(CI_95$AUC_raw_EUR_boot_Upper,CI_95$AUC_raw_SAS_boot_Upper,CI_95$AUC_raw_AFR_boot_Upper,CI_95$AUC_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower,CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper,CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper),
                    AUC_adjusted_Lower_95 = c(CI_95$AUC_adjusted_EUR_boot_Lower,CI_95$AUC_adjusted_SAS_boot_Lower,CI_95$AUC_adjusted_AFR_boot_Lower,CI_95$AUC_adjusted_AMR_boot_Lower),
                    AUC_adjusted_Upper_95 = c(CI_95$AUC_adjusted_EUR_boot_Upper,CI_95$AUC_adjusted_SAS_boot_Upper,CI_95$AUC_adjusted_AFR_boot_Upper,CI_95$AUC_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","AUC_adjusted","AUC_adjusted_Lower_95","AUC_adjusted_Upper_95")]
colnames(full_results) <- c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_AUC_adjusted","R2_AUC_adjusted_Lower_95","R2_AUC_adjusted_Upper_95")

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_Imputed <- rbind(final_results_Imputed,full_results)

write.xlsx(final_results_Imputed, file = "Metrics_Table.xlsx",sheetName = "UKB Imputed", append = TRUE,row.names = FALSE)





rm(list = ls())

final_results_WGS <- NULL

full_results <- NULL
full_results_Boot <- NULL

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  CT_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/CT/",trait,"_Bootstraps.csv"))
  CT_Boot_Results$Method <- "CT"
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  LDPred2_Results$Method <- "LDpred2"
  LDPred2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"_Bootstraps_LDPred2.csv"))
  LDPred2_Boot_Results$Method <- "LDpred2"
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  LASSOSUM2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"_Bootstraps_LASSOSum.csv"))
  LASSOSUM2_Boot_Results$Method <- "Lassosum2"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_CV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"_Bootstraps.csv"))
  RICE_CV_Boot_Results$Method <- "RICE-CV"
  colnames(RICE_CV_Boot_Results) <- colnames(CT_Boot_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  RICE_RV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"_Bootstraps.csv"))
  RICE_RV_Boot_Results$Method <- "RICE-RV"
  colnames(RICE_RV_Boot_Results) <- colnames(CT_Boot_Results)
  full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,RICE_CV_Results,RICE_RV_Results))
  full_results_Boot <- rbind(full_results_Boot,rbind(CT_Boot_Results,LDPred2_Boot_Results,LASSOSUM2_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_Boot$trait <- factor(full_results_Boot$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_SAS_boot_Lower,CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_SAS_boot_Upper,CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper),
                    R2_raw_Lower_95 = c(CI_95$R2_raw_EUR_boot_Lower,CI_95$R2_raw_SAS_boot_Lower,CI_95$R2_raw_AFR_boot_Lower,CI_95$R2_raw_AMR_boot_Lower),
                    R2_raw_Upper_95 = c(CI_95$R2_raw_EUR_boot_Upper,CI_95$R2_raw_SAS_boot_Upper,CI_95$R2_raw_AFR_boot_Upper,CI_95$R2_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower,CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper,CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper),
                    R2_adjusted_Lower_95 = c(CI_95$R2_adjusted_EUR_boot_Lower,CI_95$R2_adjusted_SAS_boot_Lower,CI_95$R2_adjusted_AFR_boot_Lower,CI_95$R2_adjusted_AMR_boot_Lower),
                    R2_adjusted_Upper_95 = c(CI_95$R2_adjusted_EUR_boot_Upper,CI_95$R2_adjusted_SAS_boot_Upper,CI_95$R2_adjusted_AFR_boot_Upper,CI_95$R2_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_adjusted","R2_adjusted_Lower_95","R2_adjusted_Upper_95")]
colnames(full_results) <- c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_AUC_adjusted","R2_AUC_adjusted_Lower_95","R2_AUC_adjusted_Upper_95")

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_WGS <- rbind(final_results_WGS,full_results)

full_results <- NULL
full_results_Boot <- NULL

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  CT_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/CT/",trait,"_Bootstraps.csv"))
  CT_Boot_Results$Method <- "CT"
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  LDPred2_Results$Method <- "LDpred2"
  LDPred2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"_Bootstraps_LDPred2.csv"))
  LDPred2_Boot_Results$Method <- "LDpred2"
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  LASSOSUM2_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"_Bootstraps_LASSOSum.csv"))
  LASSOSUM2_Boot_Results$Method <- "Lassosum2"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_CV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/CV_",trait,"_Bootstraps.csv"))
  RICE_CV_Boot_Results$Method <- "RICE-CV"
  colnames(RICE_CV_Boot_Results) <- colnames(CT_Boot_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  RICE_RV_Boot_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/RV_",trait,"_Bootstraps.csv"))
  RICE_RV_Boot_Results$Method <- "RICE-RV"
  colnames(RICE_RV_Boot_Results) <- colnames(CT_Boot_Results)
  full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,RICE_CV_Results,RICE_RV_Results))
  full_results_Boot <- rbind(full_results_Boot,rbind(CT_Boot_Results,LDPred2_Boot_Results,LASSOSUM2_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait <- factor(full_results$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_SAS_boot_Lower,CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_SAS_boot_Upper,CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper),
                    AUC_raw_Lower_95 = c(CI_95$AUC_raw_EUR_boot_Lower,CI_95$AUC_raw_SAS_boot_Lower,CI_95$AUC_raw_AFR_boot_Lower,CI_95$AUC_raw_AMR_boot_Lower),
                    AUC_raw_Upper_95 = c(CI_95$AUC_raw_EUR_boot_Upper,CI_95$AUC_raw_SAS_boot_Upper,CI_95$AUC_raw_AFR_boot_Upper,CI_95$AUC_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower,CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper,CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper),
                    AUC_adjusted_Lower_95 = c(CI_95$AUC_adjusted_EUR_boot_Lower,CI_95$AUC_adjusted_SAS_boot_Lower,CI_95$AUC_adjusted_AFR_boot_Lower,CI_95$AUC_adjusted_AMR_boot_Lower),
                    AUC_adjusted_Upper_95 = c(CI_95$AUC_adjusted_EUR_boot_Upper,CI_95$AUC_adjusted_SAS_boot_Upper,CI_95$AUC_adjusted_AFR_boot_Upper,CI_95$AUC_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","AUC_adjusted","AUC_adjusted_Lower_95","AUC_adjusted_Upper_95")]
colnames(full_results) <- c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_AUC_adjusted","R2_AUC_adjusted_Lower_95","R2_AUC_adjusted_Upper_95")

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_WGS <- rbind(final_results_WGS,full_results)

write.xlsx(final_results_WGS, file = "Metrics_Table.xlsx",sheetName = "UKB WGS", append = TRUE,row.names = FALSE)



rm(list = ls())

CTSLEB_Results <- read.csv("/data/williamsjacr/AoU_Results/CTSLEB_Results.csv")
CTSLEB_Results$Method <- "CT-SLEB"
CTSLEB_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CTSLEB_Boot.csv")
CTSLEB_Boot_Results$Method <- "CT-SLEB"
PROSPER_Results <- read.csv("/data/williamsjacr/AoU_Results/PROSPER_Results.csv")
PROSPER_Results$Method <- "PROSPER"
PROSPER_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/PROSPER_Boot.csv")
PROSPER_Boot_Results$Method <- "PROSPER"
JointPRS_Results <- read.csv("/data/williamsjacr/AoU_Results/JointPRS_Results.csv")
JointPRS_Results$Method <- "JointPRS"
JointPRS_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/JointPRS_Boot.csv")
JointPRS_Boot_Results$Method <- "JointPRS"
RICE_CV_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Results.csv")
RICE_CV_Results$Method <- "RICE-CV"
RICE_CV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Boot.csv")
RICE_CV_Boot_Results$Method <- "RICE-CV"
colnames(RICE_CV_Boot_Results) <- colnames(CTSLEB_Boot_Results)
RICE_RV_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Results.csv")
RICE_RV_Results$Method <- "RICE-RV"
RICE_RV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Boot.csv")
RICE_RV_Boot_Results$Method <- "RICE-RV"
colnames(RICE_RV_Boot_Results) <- colnames(CTSLEB_Boot_Results)

full_results <- rbind(CTSLEB_Results,PROSPER_Results,JointPRS_Results,RICE_CV_Results,RICE_RV_Results)
full_results_Boot <- rbind(CTSLEB_Boot_Results,PROSPER_Boot_Results,JointPRS_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results)

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_Boot$trait <- factor(full_results_Boot$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EAS","EUR","MID","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower,CI_95$beta_raw_EAS_boot_Lower,CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_MID_boot_Lower,CI_95$beta_raw_SAS_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper,CI_95$beta_raw_EAS_boot_Upper,CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_MID_boot_Upper,CI_95$beta_raw_SAS_boot_Upper),
                    R2_raw_Lower_95 = c(CI_95$R2_raw_AFR_boot_Lower,CI_95$R2_raw_AMR_boot_Lower,CI_95$R2_raw_EAS_boot_Lower,CI_95$R2_raw_EUR_boot_Lower,CI_95$R2_raw_MID_boot_Lower,CI_95$R2_raw_SAS_boot_Lower),
                    R2_raw_Upper_95 = c(CI_95$R2_raw_AFR_boot_Upper,CI_95$R2_raw_AMR_boot_Upper,CI_95$R2_raw_EAS_boot_Upper,CI_95$R2_raw_EUR_boot_Upper,CI_95$R2_raw_MID_boot_Upper,CI_95$R2_raw_SAS_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower,CI_95$beta_adjusted_EAS_boot_Lower,CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_MID_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper,CI_95$beta_adjusted_EAS_boot_Upper,CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_MID_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper),
                    R2_adjusted_Lower_95 = c(CI_95$R2_adjusted_AFR_boot_Lower,CI_95$R2_adjusted_AMR_boot_Lower,CI_95$R2_adjusted_EAS_boot_Lower,CI_95$R2_adjusted_EUR_boot_Lower,CI_95$R2_adjusted_MID_boot_Lower,CI_95$R2_adjusted_SAS_boot_Lower),
                    R2_adjusted_Upper_95 = c(CI_95$R2_adjusted_AFR_boot_Upper,CI_95$R2_adjusted_AMR_boot_Upper,CI_95$R2_adjusted_EAS_boot_Upper,CI_95$R2_adjusted_EUR_boot_Upper,CI_95$R2_adjusted_MID_boot_Upper,CI_95$R2_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

full_results <- full_results[,c("trait","ancestry","Method","beta_adjusted","beta_adjusted_Lower_95","beta_adjusted_Upper_95","R2_adjusted","R2_adjusted_Lower_95","R2_adjusted_Upper_95")]

full_results <- full_results[order(full_results[,1], full_results[,2]), ]

final_results_AoU <- full_results

write.xlsx(final_results_AoU, file = "Metrics_Table.xlsx",sheetName = "AoU", append = TRUE,row.names = FALSE)