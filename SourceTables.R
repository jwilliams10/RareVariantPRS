rm(list = ls())

#################################################################
### Figure 3 + Supplementary Figures 1
#################################################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)

index_mat <- NULL

causalprop_vec <- c(0.2,0.05,0.01,0.001,0.0005)
scale <- c(0,1)

count <- 1

for(j in 1:length(causalprop_vec)){
  for(q in 1:length(scale)){
    for(l in 1:100){
      index_mat <- rbind(index_mat,data.frame(i = count,Causal_Prop = causalprop_vec[j],Scale = scale[q]))
      count <- count + 1
    }
  }
}

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_70_NewProp <- NULL
results_70_NewProp_CIs <- NULL
results_70_NewProp_Comparison_CIs <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv")))
  
  Bootstraps_Comparison <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation5/Results/Common_plus_RareVariants/Comparison_Bootstraps",i,".csv")))
  
  lower_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  lower_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_95 = lower_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_95 = lower_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_95 = lower_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_95 = lower_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  upper_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_95 = upper_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_95 = upper_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_95 = upper_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_95 = upper_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  lower_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  lower_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_99 = lower_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_99 = lower_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_99 = lower_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_99 = lower_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  upper_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_99 = upper_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_99 = upper_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_99 = upper_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_99 = upper_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_70_NewProp_CIs <- rbind(results_70_NewProp_CIs,CIs_tmp)
  
  lower_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  upper_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  CI_95 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  lower_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  upper_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  CI_99 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  CIs_tmp <- inner_join(CI_95,CI_99)
  results_70_NewProp_Comparison_CIs <- rbind(results_70_NewProp_Comparison_CIs,CIs_tmp)
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_70_NewProp <- rbind(results_70_NewProp,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_70_Comparison_CIs",
                          "results_35","results_35_CIs","results_35_Comparison_CIs",
                          "results_rareprop_70","results_rareprop_70_CI","results_rareprop_70_Comparison_CIs",
                          "results_rareprop_35","results_rareprop_35_CI","results_rareprop_35_Comparison_CIs",
                          "results_70_NewProp","results_70_NewProp_CIs","results_70_NewProp_Comparison_CIs",
                          "i","Y_train","index_mat")))
}

results_70_NewProp <- inner_join(results_70_NewProp,index_mat)
results_70_NewProp$Causal_Prop <- as.character(results_70_NewProp$Causal_Prop)
results_70_NewProp$Causal_Prop[results_70_NewProp$Causal_Prop == "5e-04"] <- "0.0005"
results_70_NewProp$Causal_Prop <- paste0("Causal Prop. ",results_70_NewProp$Causal_Prop)

results_70_NewProp$Scale <- as.character(results_70_NewProp$Scale)
results_70_NewProp$Scale[results_70_NewProp$Scale == "0"] <- "Unscaled"
results_70_NewProp$Scale[results_70_NewProp$Scale == "1"] <- "Scaled"

results_70_NewProp <- data.frame(Scale = results_70_NewProp$Scale, Causal_Prop = results_70_NewProp$Causal_Prop, Method = results_70_NewProp$Method,Ancestry = results_70_NewProp$ancestry,
                                 Beta = results_70_NewProp$beta_adjusted,SE_Beta = results_70_NewProp$beta_se_adjusted,R2 = results_70_NewProp$R2_adjusted,SE_R2 = results_70_NewProp$R2_se_adjusted)
results_70_NewProp$Train_Size <- nrow(Y_train[[1]])

results_70_NewProp_CIs <- inner_join(results_70_NewProp_CIs,index_mat)
results_70_NewProp_CIs$Causal_Prop <- as.character(results_70_NewProp_CIs$Causal_Prop)
results_70_NewProp_CIs$Causal_Prop[results_70_NewProp_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_70_NewProp_CIs$Causal_Prop <- paste0("Causal Prop. ",results_70_NewProp_CIs$Causal_Prop)

results_70_NewProp_CIs$Scale <- as.character(results_70_NewProp_CIs$Scale)
results_70_NewProp_CIs$Scale[results_70_NewProp_CIs$Scale == "0"] <- "Unscaled"
results_70_NewProp_CIs$Scale[results_70_NewProp_CIs$Scale == "1"] <- "Scaled"

results_70_NewProp_CIs$Train_Size <- nrow(Y_train[[1]])

results_70_NewProp_Comparison_CIs <- inner_join(results_70_NewProp_Comparison_CIs,index_mat)
results_70_NewProp_Comparison_CIs$Causal_Prop <- as.character(results_70_NewProp_Comparison_CIs$Causal_Prop)
results_70_NewProp_Comparison_CIs$Causal_Prop[results_70_NewProp_Comparison_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_70_NewProp_Comparison_CIs$Causal_Prop <- paste0("Causal Prop. ",results_70_NewProp_Comparison_CIs$Causal_Prop)

results_70_NewProp_Comparison_CIs$Scale <- as.character(results_70_NewProp_Comparison_CIs$Scale)
results_70_NewProp_Comparison_CIs$Scale[results_70_NewProp_Comparison_CIs$Scale == "0"] <- "Unscaled"
results_70_NewProp_Comparison_CIs$Scale[results_70_NewProp_Comparison_CIs$Scale == "1"] <- "Scaled"

results_70_NewProp_Comparison_CIs$Train_Size <- nrow(Y_train[[1]])



load("/data/williamsjacr/UKB_WES_Simulation/Simulation6/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_35_NewProp <- NULL
results_35_NewProp_CIs <- NULL
results_35_NewProp_Comparison_CIs <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv")))
  
  Bootstraps_Comparison <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation6/Results/Common_plus_RareVariants/Comparison_Bootstraps",i,".csv")))
  
  lower_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  lower_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_95 = lower_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_95 = lower_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_95 = lower_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_95 = lower_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  upper_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_95 = upper_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_95 = upper_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_95 = upper_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_95 = upper_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  lower_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  lower_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_99 = lower_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_99 = lower_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_99 = lower_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_99 = lower_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  upper_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_99 = upper_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_99 = upper_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_99 = upper_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_99 = upper_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_35_NewProp_CIs <- rbind(results_35_NewProp_CIs,CIs_tmp)
  
  lower_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  upper_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  CI_95 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  lower_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  upper_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  CI_99 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  CIs_tmp <- inner_join(CI_95,CI_99)
  results_35_NewProp_Comparison_CIs <- rbind(results_35_NewProp_Comparison_CIs,CIs_tmp)
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_35_NewProp <- rbind(results_35_NewProp,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_70_Comparison_CIs",
                          "results_35","results_35_CIs","results_35_Comparison_CIs",
                          "results_rareprop_70","results_rareprop_70_CI","results_rareprop_70_Comparison_CIs",
                          "results_rareprop_35","results_rareprop_35_CI","results_rareprop_35_Comparison_CIs",
                          "results_70_NewProp","results_70_NewProp_CIs","results_70_NewProp_Comparison_CIs",
                          "results_35_NewProp","results_35_NewProp_CIs","results_35_NewProp_Comparison_CIs",
                          "i","Y_train","index_mat")))
}

results_35_NewProp <- inner_join(results_35_NewProp,index_mat)
results_35_NewProp$Causal_Prop <- as.character(results_35_NewProp$Causal_Prop)
results_35_NewProp$Causal_Prop[results_35_NewProp$Causal_Prop == "5e-04"] <- "0.0005"
results_35_NewProp$Causal_Prop <- paste0("Causal Prop. ",results_35_NewProp$Causal_Prop)

results_35_NewProp$Scale <- as.character(results_35_NewProp$Scale)
results_35_NewProp$Scale[results_35_NewProp$Scale == "0"] <- "Unscaled"
results_35_NewProp$Scale[results_35_NewProp$Scale == "1"] <- "Scaled"

results_35_NewProp <- data.frame(Scale = results_35_NewProp$Scale, Causal_Prop = results_35_NewProp$Causal_Prop, Method = results_35_NewProp$Method,Ancestry = results_35_NewProp$ancestry,
                                 Beta = results_35_NewProp$beta_adjusted,SE_Beta = results_35_NewProp$beta_se_adjusted,R2 = results_35_NewProp$R2_adjusted,SE_R2 = results_35_NewProp$R2_se_adjusted)
results_35_NewProp$Train_Size <- nrow(Y_train[[1]])

results_35_NewProp_CIs <- inner_join(results_35_NewProp_CIs,index_mat)
results_35_NewProp_CIs$Causal_Prop <- as.character(results_35_NewProp_CIs$Causal_Prop)
results_35_NewProp_CIs$Causal_Prop[results_35_NewProp_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_35_NewProp_CIs$Causal_Prop <- paste0("Causal Prop. ",results_35_NewProp_CIs$Causal_Prop)

results_35_NewProp_CIs$Scale <- as.character(results_35_NewProp_CIs$Scale)
results_35_NewProp_CIs$Scale[results_35_NewProp_CIs$Scale == "0"] <- "Unscaled"
results_35_NewProp_CIs$Scale[results_35_NewProp_CIs$Scale == "1"] <- "Scaled"

results_35_NewProp_CIs$Train_Size <- nrow(Y_train[[1]])

results_35_NewProp_Comparison_CIs <- inner_join(results_35_NewProp_Comparison_CIs,index_mat)
results_35_NewProp_Comparison_CIs$Causal_Prop <- as.character(results_35_NewProp_Comparison_CIs$Causal_Prop)
results_35_NewProp_Comparison_CIs$Causal_Prop[results_35_NewProp_Comparison_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_35_NewProp_Comparison_CIs$Causal_Prop <- paste0("Causal Prop. ",results_35_NewProp_Comparison_CIs$Causal_Prop)

results_35_NewProp_Comparison_CIs$Scale <- as.character(results_35_NewProp_Comparison_CIs$Scale)
results_35_NewProp_Comparison_CIs$Scale[results_35_NewProp_Comparison_CIs$Scale == "0"] <- "Unscaled"
results_35_NewProp_Comparison_CIs$Scale[results_35_NewProp_Comparison_CIs$Scale == "1"] <- "Scaled"

results_35_NewProp_Comparison_CIs$Train_Size <- nrow(Y_train[[1]])



load("/data/williamsjacr/UKB_WES_Simulation/Simulation7/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_70_NewProp_RareProp <- NULL
results_70_NewProp_RareProp_CIs <- NULL
results_70_NewProp_RareProp_Comparison_CIs <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv")))
  
  Bootstraps_Comparison <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation7/Results/Common_plus_RareVariants/Comparison_Bootstraps",i,".csv")))
  
  lower_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  lower_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_95 = lower_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_95 = lower_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_95 = lower_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_95 = lower_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  upper_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_95 = upper_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_95 = upper_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_95 = upper_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_95 = upper_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  lower_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  lower_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_99 = lower_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_99 = lower_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_99 = lower_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_99 = lower_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  upper_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_99 = upper_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_99 = upper_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_99 = upper_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_99 = upper_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_70_NewProp_RareProp_CIs <- rbind(results_70_NewProp_RareProp_CIs,CIs_tmp)
  
  lower_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  upper_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  CI_95 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  lower_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  upper_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  CI_99 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  CIs_tmp <- inner_join(CI_95,CI_99)
  results_70_NewProp_RareProp_Comparison_CIs <- rbind(results_70_NewProp_RareProp_Comparison_CIs,CIs_tmp)
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_70_NewProp_RareProp <- rbind(results_70_NewProp_RareProp,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_70_Comparison_CIs",
                          "results_35","results_35_CIs","results_35_Comparison_CIs",
                          "results_rareprop_70","results_rareprop_70_CI","results_rareprop_70_Comparison_CIs",
                          "results_rareprop_35","results_rareprop_35_CI","results_rareprop_35_Comparison_CIs",
                          "results_70_NewProp","results_70_NewProp_CIs","results_70_NewProp_Comparison_CIs",
                          "results_35_NewProp","results_35_NewProp_CIs","results_35_NewProp_Comparison_CIs",
                          "results_70_NewProp_RareProp","results_70_NewProp_RareProp_CIs","results_70_NewProp_RareProp_Comparison_CIs",
                          "i","Y_train","index_mat")))
}

results_70_NewProp_RareProp <- inner_join(results_70_NewProp_RareProp,index_mat)
results_70_NewProp_RareProp$Causal_Prop <- as.character(results_70_NewProp_RareProp$Causal_Prop)
results_70_NewProp_RareProp$Causal_Prop[results_70_NewProp_RareProp$Causal_Prop == "5e-04"] <- "0.0005"
results_70_NewProp_RareProp$Causal_Prop <- paste0("Causal Prop. ",results_70_NewProp_RareProp$Causal_Prop)

results_70_NewProp_RareProp$Scale <- as.character(results_70_NewProp_RareProp$Scale)
results_70_NewProp_RareProp$Scale[results_70_NewProp_RareProp$Scale == "0"] <- "Unscaled"
results_70_NewProp_RareProp$Scale[results_70_NewProp_RareProp$Scale == "1"] <- "Scaled"

results_70_NewProp_RareProp <- data.frame(Scale = results_70_NewProp_RareProp$Scale, Causal_Prop = results_70_NewProp_RareProp$Causal_Prop, Method = results_70_NewProp_RareProp$Method,Ancestry = results_70_NewProp_RareProp$ancestry,
                                          Beta = results_70_NewProp_RareProp$beta_adjusted,SE_Beta = results_70_NewProp_RareProp$beta_se_adjusted,R2 = results_70_NewProp_RareProp$R2_adjusted,SE_R2 = results_70_NewProp_RareProp$R2_se_adjusted)
results_70_NewProp_RareProp$Train_Size <- nrow(Y_train[[1]])

results_70_NewProp_RareProp_CIs <- inner_join(results_70_NewProp_RareProp_CIs,index_mat)
results_70_NewProp_RareProp_CIs$Causal_Prop <- as.character(results_70_NewProp_RareProp_CIs$Causal_Prop)
results_70_NewProp_RareProp_CIs$Causal_Prop[results_70_NewProp_RareProp_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_70_NewProp_RareProp_CIs$Causal_Prop <- paste0("Causal Prop. ",results_70_NewProp_RareProp_CIs$Causal_Prop)

results_70_NewProp_RareProp_CIs$Scale <- as.character(results_70_NewProp_RareProp_CIs$Scale)
results_70_NewProp_RareProp_CIs$Scale[results_70_NewProp_RareProp_CIs$Scale == "0"] <- "Unscaled"
results_70_NewProp_RareProp_CIs$Scale[results_70_NewProp_RareProp_CIs$Scale == "1"] <- "Scaled"

results_70_NewProp_RareProp_CIs$Train_Size <- nrow(Y_train[[1]])

results_70_NewProp_RareProp_Comparison_CIs <- inner_join(results_70_NewProp_RareProp_Comparison_CIs,index_mat)
results_70_NewProp_RareProp_Comparison_CIs$Causal_Prop <- as.character(results_70_NewProp_RareProp_Comparison_CIs$Causal_Prop)
results_70_NewProp_RareProp_Comparison_CIs$Causal_Prop[results_70_NewProp_RareProp_Comparison_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_70_NewProp_RareProp_Comparison_CIs$Causal_Prop <- paste0("Causal Prop. ",results_70_NewProp_RareProp_Comparison_CIs$Causal_Prop)

results_70_NewProp_RareProp_Comparison_CIs$Scale <- as.character(results_70_NewProp_RareProp_Comparison_CIs$Scale)
results_70_NewProp_RareProp_Comparison_CIs$Scale[results_70_NewProp_RareProp_Comparison_CIs$Scale == "0"] <- "Unscaled"
results_70_NewProp_RareProp_Comparison_CIs$Scale[results_70_NewProp_RareProp_Comparison_CIs$Scale == "1"] <- "Scaled"

results_70_NewProp_RareProp_Comparison_CIs$Train_Size <- nrow(Y_train[[1]])



load("/data/williamsjacr/UKB_WES_Simulation/Simulation8/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_35_NewProp_RareProp <- NULL
results_35_NewProp_RareProp_CIs <- NULL
results_35_NewProp_RareProp_Comparison_CIs <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv")))
  
  Bootstraps_Comparison <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation8/Results/Common_plus_RareVariants/Comparison_Bootstraps",i,".csv")))
  
  lower_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  lower_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_95 = lower_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_95 = lower_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_95 = lower_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_95 = lower_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  upper_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_95 = upper_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_95 = upper_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_95 = upper_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_95 = upper_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  lower_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  lower_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_99 = lower_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_99 = lower_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_99 = lower_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_99 = lower_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  upper_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry =c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_99 = upper_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_99 = upper_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_99 = upper_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_99 = upper_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_35_NewProp_RareProp_CIs <- rbind(results_35_NewProp_RareProp_CIs,CIs_tmp)
  
  lower_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  upper_95 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  CI_95 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_95 = lower_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_95 = upper_95[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  lower_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  upper_99 <- apply(Bootstraps_Comparison,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  CI_99 <- data.frame(i = Bootstraps_Comparison$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),
                      R2_raw_RICE_vs_CT_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_CT_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_CT","R2_raw_SAS_RICE_vs_CT","R2_raw_AMR_RICE_vs_CT","R2_raw_AFR_RICE_vs_CT")],
                      R2_raw_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_LDpred2","R2_raw_SAS_RICE_vs_LDpred2","R2_raw_AMR_RICE_vs_LDpred2","R2_raw_AFR_RICE_vs_LDpred2")],
                      R2_raw_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_raw_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_raw_EUR_RICE_vs_Lassosum2","R2_raw_SAS_RICE_vs_Lassosum2","R2_raw_AMR_RICE_vs_Lassosum2","R2_raw_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_CT_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_CT_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_CT","R2_adjusted_SAS_RICE_vs_CT","R2_adjusted_AMR_RICE_vs_CT","R2_adjusted_AFR_RICE_vs_CT")],
                      R2_adjusted_RICE_vs_LDpred2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_LDpred2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_LDpred2","R2_adjusted_SAS_RICE_vs_LDpred2","R2_adjusted_AMR_RICE_vs_LDpred2","R2_adjusted_AFR_RICE_vs_LDpred2")],
                      R2_adjusted_RICE_vs_Lassosum2_Lower_99 = lower_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")],
                      R2_adjusted_RICE_vs_Lassosum2_Upper_99 = upper_99[c("R2_adjusted_EUR_RICE_vs_Lassosum2","R2_adjusted_SAS_RICE_vs_Lassosum2","R2_adjusted_AMR_RICE_vs_Lassosum2","R2_adjusted_AFR_RICE_vs_Lassosum2")])
  
  CIs_tmp <- inner_join(CI_95,CI_99)
  results_35_NewProp_RareProp_Comparison_CIs <- rbind(results_35_NewProp_RareProp_Comparison_CIs,CIs_tmp)
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_35_NewProp_RareProp <- rbind(results_35_NewProp_RareProp,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_70_Comparison_CIs",
                          "results_35","results_35_CIs","results_35_Comparison_CIs",
                          "results_rareprop_70","results_rareprop_70_CI","results_rareprop_70_Comparison_CIs",
                          "results_rareprop_35","results_rareprop_35_CI","results_rareprop_35_Comparison_CIs",
                          "results_70_NewProp","results_70_NewProp_CIs","results_70_NewProp_Comparison_CIs",
                          "results_35_NewProp","results_35_NewProp_CIs","results_35_NewProp_Comparison_CIs",
                          "results_70_NewProp_RareProp","results_70_NewProp_RareProp_CIs","results_70_NewProp_RareProp_Comparison_CIs",
                          "results_35_NewProp_RareProp","results_35_NewProp_RareProp_CIs","results_35_NewProp_RareProp_Comparison_CIs",
                          "i","Y_train","index_mat")))
}

results_35_NewProp_RareProp <- inner_join(results_35_NewProp_RareProp,index_mat)
results_35_NewProp_RareProp$Causal_Prop <- as.character(results_35_NewProp_RareProp$Causal_Prop)
results_35_NewProp_RareProp$Causal_Prop[results_35_NewProp_RareProp$Causal_Prop == "5e-04"] <- "0.0005"
results_35_NewProp_RareProp$Causal_Prop <- paste0("Causal Prop. ",results_35_NewProp_RareProp$Causal_Prop)

results_35_NewProp_RareProp$Scale <- as.character(results_35_NewProp_RareProp$Scale)
results_35_NewProp_RareProp$Scale[results_35_NewProp_RareProp$Scale == "0"] <- "Unscaled"
results_35_NewProp_RareProp$Scale[results_35_NewProp_RareProp$Scale == "1"] <- "Scaled"

results_35_NewProp_RareProp <- data.frame(Scale = results_35_NewProp_RareProp$Scale, Causal_Prop = results_35_NewProp_RareProp$Causal_Prop, Method = results_35_NewProp_RareProp$Method,Ancestry = results_35_NewProp_RareProp$ancestry,
                                          Beta = results_35_NewProp_RareProp$beta_adjusted,SE_Beta = results_35_NewProp_RareProp$beta_se_adjusted,R2 = results_35_NewProp_RareProp$R2_adjusted,SE_R2 = results_35_NewProp_RareProp$R2_se_adjusted)
results_35_NewProp_RareProp$Train_Size <- nrow(Y_train[[1]])

results_35_NewProp_RareProp_CIs <- inner_join(results_35_NewProp_RareProp_CIs,index_mat)
results_35_NewProp_RareProp_CIs$Causal_Prop <- as.character(results_35_NewProp_RareProp_CIs$Causal_Prop)
results_35_NewProp_RareProp_CIs$Causal_Prop[results_35_NewProp_RareProp_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_35_NewProp_RareProp_CIs$Causal_Prop <- paste0("Causal Prop. ",results_35_NewProp_RareProp_CIs$Causal_Prop)

results_35_NewProp_RareProp_CIs$Scale <- as.character(results_35_NewProp_RareProp_CIs$Scale)
results_35_NewProp_RareProp_CIs$Scale[results_35_NewProp_RareProp_CIs$Scale == "0"] <- "Unscaled"
results_35_NewProp_RareProp_CIs$Scale[results_35_NewProp_RareProp_CIs$Scale == "1"] <- "Scaled"

results_35_NewProp_RareProp_CIs$Train_Size <- nrow(Y_train[[1]])

results_35_NewProp_RareProp_Comparison_CIs <- inner_join(results_35_NewProp_RareProp_Comparison_CIs,index_mat)
results_35_NewProp_RareProp_Comparison_CIs$Causal_Prop <- as.character(results_35_NewProp_RareProp_Comparison_CIs$Causal_Prop)
results_35_NewProp_RareProp_Comparison_CIs$Causal_Prop[results_35_NewProp_RareProp_Comparison_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_35_NewProp_RareProp_Comparison_CIs$Causal_Prop <- paste0("Causal Prop. ",results_35_NewProp_RareProp_Comparison_CIs$Causal_Prop)

results_35_NewProp_RareProp_Comparison_CIs$Scale <- as.character(results_35_NewProp_RareProp_Comparison_CIs$Scale)
results_35_NewProp_RareProp_Comparison_CIs$Scale[results_35_NewProp_RareProp_Comparison_CIs$Scale == "0"] <- "Unscaled"
results_35_NewProp_RareProp_Comparison_CIs$Scale[results_35_NewProp_RareProp_Comparison_CIs$Scale == "1"] <- "Scaled"

results_35_NewProp_RareProp_Comparison_CIs$Train_Size <- nrow(Y_train[[1]])





results_NewProp <- rbind(results_35_NewProp,results_70_NewProp)
results_NewProp_CI <- rbind(results_35_NewProp_CIs,results_70_NewProp_CIs)
results_NewProp_Comparison_CI <- rbind(results_35_NewProp_Comparison_CIs,results_70_NewProp_Comparison_CIs)
results_NewProp_RareProp <- rbind(results_35_NewProp_RareProp,results_70_NewProp_RareProp)
results_NewProp_RareProp_CI <- rbind(results_35_NewProp_RareProp_CIs,results_70_NewProp_RareProp_CIs)
results_NewProp_RareProp_Comparison_CI <- rbind(results_35_NewProp_RareProp_Comparison_CIs,results_70_NewProp_RareProp_Comparison_CIs)

results_NewProp$Train_Size <- format(results_NewProp$Train_Size,big.mark=",", trim=TRUE)
results_NewProp_RareProp$Train_Size <- format(results_NewProp_RareProp$Train_Size,big.mark=",", trim=TRUE)

results_NewProp_CI$Train_Size <- format(results_NewProp_CI$Train_Size,big.mark=",", trim=TRUE)
results_NewProp_RareProp_CI$Train_Size <- format(results_NewProp_RareProp_CI$Train_Size,big.mark=",", trim=TRUE)

results_NewProp_Comparison_CI$Train_Size <- format(results_NewProp_Comparison_CI$Train_Size,big.mark=",", trim=TRUE)
results_NewProp_RareProp_Comparison_CI$Train_Size <- format(results_NewProp_RareProp_Comparison_CI$Train_Size,big.mark=",", trim=TRUE)

results_NewProp$Train_Size <- paste0("n = ",results_NewProp$Train_Size)
results_NewProp_RareProp$Train_Size <- paste0("n = ",results_NewProp_RareProp$Train_Size)

results_NewProp_CI$Train_Size <- paste0("n = ",results_NewProp_CI$Train_Size)
results_NewProp_RareProp_CI$Train_Size <- paste0("n = ",results_NewProp_RareProp_CI$Train_Size)

results_NewProp_Comparison_CI$Train_Size <- paste0("n = ",results_NewProp_Comparison_CI$Train_Size)
results_NewProp_RareProp_Comparison_CI$Train_Size <- paste0("n = ",results_NewProp_RareProp_Comparison_CI$Train_Size)

rm(list=setdiff(ls(), c("results","results_CI","results_rareprop","results_rareprop_CI","results_NewProp","results_NewProp_CI","results_Comparisons_CI","results_rareprop_Comparison_CI","results_NewProp_Comparison_CI","results_NewProp_RareProp","results_NewProp_RareProp_CI","results_NewProp_RareProp_Comparison_CI")))


results_NewProp$Beta[results_NewProp$Method %in% c("LDPred","LASSOSum")] <- -1*results_NewProp$Beta[results_NewProp$Method %in% c("LDPred","LASSOSum")]
results_NewProp$Beta[results_NewProp$Beta < 0] <- 0
results_NewProp$Beta[results_NewProp$Beta > 1] <- 0

results_NewProp_RareProp$Beta[results_NewProp_RareProp$Method %in% c("LDPred","LASSOSum")] <- -1*results_NewProp_RareProp$Beta[results_NewProp_RareProp$Method %in% c("LDPred","LASSOSum")]
results_NewProp_RareProp$Beta[results_NewProp_RareProp$Beta < 0] <- 0
results_NewProp_RareProp$Beta[results_NewProp_RareProp$Beta > 1] <- 0

results_NewProp <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results_NewProp,mean)
results_NewProp_RareProp <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results_NewProp_RareProp,mean)

results_NewProp_CI <- aggregate(.~Scale + Causal_Prop + Train_Size + Ancestry,data = subset(results_NewProp_CI,select = -c(i)),mean)
results_NewProp_RareProp_CI <- aggregate(.~Scale + Causal_Prop + Train_Size + Ancestry,data = subset(results_NewProp_RareProp_CI,select = -c(i)),mean)

results_NewProp_Comparison_CI <- aggregate(.~Scale + Causal_Prop + Train_Size + Ancestry,data = subset(results_NewProp_Comparison_CI,select = -c(i)),mean)
results_NewProp_RareProp_Comparison_CI <- aggregate(.~Scale + Causal_Prop + Train_Size + Ancestry,data = subset(results_NewProp_RareProp_Comparison_CI,select = -c(i)),mean)

overall_results_NewProp <- results_NewProp[results_NewProp$Method %in% c("CT","LASSOSum","LDPred","RICE-CV","RICE-RV"),]
overall_results_NewProp_RareProp <- results_NewProp_RareProp[results_NewProp_RareProp$Method %in% c("CT","LASSOSum","LDPred","RICE-CV","RICE-RV"),]

overall_results_NewProp <- overall_results_NewProp[overall_results_NewProp$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results_NewProp$Method[overall_results_NewProp$Method == "RICE-CV"] <- "RICE-CV" 
overall_results_NewProp$Method[overall_results_NewProp$Method == "RICE-RV"] <- "RICE-RV" 
overall_results_NewProp$Method[overall_results_NewProp$Method == "LDPred"] <- "LDpred2"
overall_results_NewProp$Method[overall_results_NewProp$Method == "LASSOSum"] <- "Lassosum2"

overall_results_NewProp$Method1 <- overall_results_NewProp$Method
overall_results_NewProp$Method <- factor(overall_results_NewProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
overall_results_NewProp$Method1[overall_results_NewProp$Method1 == "RICE-RV"] <- "RICE-CV"
overall_results_NewProp$Method1 <- factor(overall_results_NewProp$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))
overall_results_NewProp$Ancestry <- factor(overall_results_NewProp$Ancestry,levels = c("AFR","AMR","EUR","SAS"))

overall_results_NewProp_RareProp <- overall_results_NewProp_RareProp[overall_results_NewProp_RareProp$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results_NewProp_RareProp$Method[overall_results_NewProp_RareProp$Method == "RICE-CV"] <- "RICE-CV" 
overall_results_NewProp_RareProp$Method[overall_results_NewProp_RareProp$Method == "RICE-RV"] <- "RICE-RV" 
overall_results_NewProp_RareProp$Method[overall_results_NewProp_RareProp$Method == "LDPred"] <- "LDpred2"
overall_results_NewProp_RareProp$Method[overall_results_NewProp_RareProp$Method == "LASSOSum"] <- "Lassosum2"

overall_results_NewProp_RareProp$Method1 <- overall_results_NewProp_RareProp$Method
overall_results_NewProp_RareProp$Method <- factor(overall_results_NewProp_RareProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
overall_results_NewProp_RareProp$Method1[overall_results_NewProp_RareProp$Method1 == "RICE-RV"] <- "RICE-CV"
overall_results_NewProp_RareProp$Method1 <- factor(overall_results_NewProp_RareProp$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))
overall_results_NewProp_RareProp$Ancestry <- factor(overall_results_NewProp_RareProp$Ancestry,levels = c("AFR","AMR","EUR","SAS"))

results_NewProp_CI$Method <- "RICE-RV"
overall_results_NewProp <- left_join(overall_results_NewProp,results_NewProp_CI)
overall_results_NewProp$Method <- factor(overall_results_NewProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

overall_results_NewProp$group1 <- "RICE-CV"
overall_results_NewProp$group2 <- "RICE-CV"
overall_results_NewProp$p.signif_beta <- ""
overall_results_NewProp$p.signif_beta[overall_results_NewProp$Method == "RICE-CV"] <- ifelse(overall_results_NewProp$Beta_Adjusted_Lower_99[overall_results_NewProp$Method == "RICE-RV"] > 0,"***",ifelse(overall_results_NewProp$Beta_Adjusted_Lower_95[overall_results_NewProp$Method == "RICE-RV"] > 0,"**",""))

overall_results_NewProp$position[overall_results_NewProp$Method == "RICE-CV"] <- overall_results_NewProp$Beta[overall_results_NewProp$Method == "RICE-CV"] + overall_results_NewProp$Beta[overall_results_NewProp$Method == "RICE-RV"] + 0.03
ylim_NewProp <- max(c(overall_results_NewProp$Beta[overall_results_NewProp$Method == "RICE-CV"] + overall_results_NewProp$Beta[overall_results_NewProp$Method == "RICE-RV"],overall_results_NewProp$Beta)) + 0.05


results_NewProp_RareProp_CI$Method <- "RICE-RV"
overall_results_NewProp_RareProp <- left_join(overall_results_NewProp_RareProp,results_NewProp_RareProp_CI)
overall_results_NewProp_RareProp$Method <- factor(overall_results_NewProp_RareProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

overall_results_NewProp_RareProp$group1 <- "RICE-CV"
overall_results_NewProp_RareProp$group2 <- "RICE-CV"
overall_results_NewProp_RareProp$p.signif_beta <- ""
overall_results_NewProp_RareProp$p.signif_beta[overall_results_NewProp_RareProp$Method == "RICE-CV"] <- ifelse(overall_results_NewProp_RareProp$Beta_Adjusted_Lower_99[overall_results_NewProp_RareProp$Method == "RICE-RV"] > 0,"***",ifelse(overall_results_NewProp_RareProp$Beta_Adjusted_Lower_95[overall_results_NewProp_RareProp$Method == "RICE-RV"] > 0,"**",""))

overall_results_NewProp_RareProp$position[overall_results_NewProp_RareProp$Method == "RICE-CV"] <- overall_results_NewProp_RareProp$Beta[overall_results_NewProp_RareProp$Method == "RICE-CV"] + overall_results_NewProp_RareProp$Beta[overall_results_NewProp_RareProp$Method == "RICE-RV"] + 0.03
ylim_NewProp_RareProp <- max(c(overall_results_NewProp_RareProp$Beta[overall_results_NewProp_RareProp$Method == "RICE-CV"] + overall_results_NewProp_RareProp$Beta[overall_results_NewProp_RareProp$Method == "RICE-RV"],overall_results_NewProp_RareProp$Beta)) + 0.05


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
            axis.title.x = element_blank(),
            axis.text.x = element_blank(), 
            axis.text.y = element_text(size = 10), 
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
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#E23838","#973999","#009cdf")), ...)
  
}

overall_results_NewProp <- overall_results_NewProp[overall_results_NewProp$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01"),]
overall_results_NewProp_RareProp <- overall_results_NewProp_RareProp[overall_results_NewProp_RareProp$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01"),]

Fig3_Betas <- overall_results_NewProp[(overall_results_NewProp$Scale == "Unscaled") & (overall_results_NewProp$Train_Size == "n = 98,343"),]
FigS1_All_Betas <- overall_results_NewProp[!((overall_results_NewProp$Scale == "Unscaled") & (overall_results_NewProp$Train_Size == "n = 98,343")),]
FigS1_Prop_Betas <- overall_results_NewProp_RareProp

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#973999","#009cdf")), ...)
}


overall_results_NewProp <- overall_results_NewProp[,c("Method","Scale","Causal_Prop","Train_Size","Ancestry","Beta","SE_Beta","R2","SE_R2")]

results_NewProp_Comparison_CI$Method <- "RICE-CV"
overall_results_NewProp <- left_join(overall_results_NewProp, results_NewProp_Comparison_CI)
overall_results_NewProp$Method <- factor(overall_results_NewProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

overall_results_NewProp$Method <- as.character(overall_results_NewProp$Method)
overall_results_NewProp <- overall_results_NewProp[overall_results_NewProp$Method != "RICE-RV",]
overall_results_NewProp$Method[overall_results_NewProp$Method == "RICE-CV"] <- "RICE"
overall_results_NewProp$Method <- factor(overall_results_NewProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

overall_results_NewProp_RareProp <- overall_results_NewProp_RareProp[,c("Method","Scale","Causal_Prop","Train_Size","Ancestry","Beta","SE_Beta","R2","SE_R2")]

results_NewProp_RareProp_Comparison_CI$Method <- "RICE-CV"
overall_results_NewProp_RareProp <- left_join(overall_results_NewProp_RareProp, results_NewProp_RareProp_Comparison_CI)
overall_results_NewProp_RareProp$Method <- factor(overall_results_NewProp_RareProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

overall_results_NewProp_RareProp$Method <- as.character(overall_results_NewProp_RareProp$Method)
overall_results_NewProp_RareProp <- overall_results_NewProp_RareProp[overall_results_NewProp_RareProp$Method != "RICE-RV",]
overall_results_NewProp_RareProp$Method[overall_results_NewProp_RareProp$Method == "RICE-CV"] <- "RICE"
overall_results_NewProp_RareProp$Method <- factor(overall_results_NewProp_RareProp$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

overall_results_NewProp$group1 <- "RICE"
overall_results_NewProp$group2 <- "RICE"
overall_results_NewProp$p.signif_beta1 <- ""
overall_results_NewProp$p.signif_beta2 <- ""

overall_results_NewProp_RareProp$group1 <- "RICE"
overall_results_NewProp_RareProp$group2 <- "RICE"
overall_results_NewProp_RareProp$p.signif_beta1 <- ""
overall_results_NewProp_RareProp$p.signif_beta2 <- ""

for(scale in c("Scaled","Unscaled")){
  for(causal_prop in c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")){
    for(train_size in c("n = 98,343","n = 49,173")){
      for(anc in c("AFR","EUR","SAS","AMR")){
        tmp <- overall_results_NewProp[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc,]
        max_R2_notRICE <- max(tmp$R2[tmp$Method != "RICE"])
        Best_Method <- tmp$Method[tmp$R2 == max_R2_notRICE]
        Improved_R2 <- round((tmp$R2[tmp$Method == "RICE"]/max_R2_notRICE - 1)*100,digits = 2)
        
        if(Best_Method == "CT"){
          overall_results_NewProp$p.signif_beta1[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc & overall_results_NewProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
          overall_results_NewProp$p.signif_beta2[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc & overall_results_NewProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0("(***)"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0("(**)"),""))
        }else if(Best_Method == "LDpred2"){
          overall_results_NewProp$p.signif_beta1[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc & overall_results_NewProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
          overall_results_NewProp$p.signif_beta2[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc & overall_results_NewProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
        }else{
          overall_results_NewProp$p.signif_beta1[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc & overall_results_NewProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
          overall_results_NewProp$p.signif_beta2[overall_results_NewProp$Scale == scale & overall_results_NewProp$Causal_Prop == causal_prop & overall_results_NewProp$Train_Size == train_size & overall_results_NewProp$Ancestry == anc & overall_results_NewProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
        }
        
        tmp <- overall_results_NewProp_RareProp[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc,]
        max_R2_notRICE <- max(tmp$R2[tmp$Method != "RICE"])
        Best_Method <- tmp$Method[tmp$R2 == max_R2_notRICE]
        Improved_R2 <- round((tmp$R2[tmp$Method == "RICE"]/max_R2_notRICE - 1)*100,digits = 2)
        
        if(Best_Method == "CT"){
          overall_results_NewProp_RareProp$p.signif_beta1[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc & overall_results_NewProp_RareProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
          overall_results_NewProp_RareProp$p.signif_beta2[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc & overall_results_NewProp_RareProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0("(***)"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0("(**)"),""))
        }else if(Best_Method == "LDpred2"){
          overall_results_NewProp_RareProp$p.signif_beta1[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc & overall_results_NewProp_RareProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
          overall_results_NewProp_RareProp$p.signif_beta2[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc & overall_results_NewProp_RareProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
        }else{
          overall_results_NewProp_RareProp$p.signif_beta1[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc & overall_results_NewProp_RareProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
          overall_results_NewProp_RareProp$p.signif_beta2[overall_results_NewProp_RareProp$Scale == scale & overall_results_NewProp_RareProp$Causal_Prop == causal_prop & overall_results_NewProp_RareProp$Train_Size == train_size & overall_results_NewProp_RareProp$Ancestry == anc & overall_results_NewProp_RareProp$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
        }
      }
    }
  }
}

overall_results_NewProp$position1 <- NA
overall_results_NewProp$position2 <- NA
overall_results_NewProp$position1[overall_results_NewProp$Method == "RICE"] <- overall_results_NewProp$R2[overall_results_NewProp$Method == "RICE"] + 0.005
overall_results_NewProp$position2[overall_results_NewProp$Method == "RICE"] <- overall_results_NewProp$R2[overall_results_NewProp$Method == "RICE"] + 0.02
ylim_NewProp <- max(c(overall_results_NewProp$R2)) + 0.03

overall_results_NewProp_RareProp$position1 <- NA
overall_results_NewProp_RareProp$position2 <- NA
overall_results_NewProp_RareProp$position1[overall_results_NewProp_RareProp$Method == "RICE"] <- overall_results_NewProp_RareProp$R2[overall_results_NewProp_RareProp$Method == "RICE"] + 0.005
overall_results_NewProp_RareProp$position2[overall_results_NewProp_RareProp$Method == "RICE"] <- overall_results_NewProp_RareProp$R2[overall_results_NewProp_RareProp$Method == "RICE"] + 0.02
ylim_NewProp_RareProp <- max(c(overall_results_NewProp_RareProp$R2)) + 0.03

Fig3_R2 <- overall_results_NewProp[(overall_results_NewProp$Scale == "Unscaled") & (overall_results_NewProp$Train_Size == "n = 98,343"),]
FigS1_All_R2 <- overall_results_NewProp[!((overall_results_NewProp$Scale == "Unscaled") & (overall_results_NewProp$Train_Size == "n = 98,343")),]
FigS1_Prop_R2 <- overall_results_NewProp_RareProp

#################################################################
### Supplementary Figure 2
#################################################################

h2_dat <- read.csv("/data/williamsjacr/UKB_WES_Simulation/Sim_Characteristics_NewProps.csv")

h2_dat_compressed <- aggregate(.~Ancestry + causal_prop + scaled,data = h2_dat[,c("Ancestry","causal_prop","scaled","h2_rare","Average_Burden_rare")],mean)
colnames(h2_dat_compressed) <- c("Ancestry","Causal_Prop","Scaled","Average_h2","Average_Burden")
h2_dat_compressed_se <- aggregate(h2_rare~Ancestry + causal_prop + scaled,data = h2_dat,function(x){quantile(x,0.025)})
colnames(h2_dat_compressed_se) <- c("Ancestry","Causal_Prop","Scaled","Q_025")
h2_dat_compressed <- inner_join(h2_dat_compressed,h2_dat_compressed_se)
h2_dat_compressed_se <- aggregate(h2_rare~Ancestry + causal_prop + scaled,data = h2_dat,function(x){quantile(x,0.975)})
colnames(h2_dat_compressed_se) <- c("Ancestry","Causal_Prop","Scaled","Q_975")
h2_dat_compressed <- inner_join(h2_dat_compressed,h2_dat_compressed_se)


h2_dat_compressed$Causal_Prop <- paste0("Causal Prop. ",h2_dat_compressed$Causal_Prop)
h2_dat_compressed$Scaled <- paste0("Scaled: ",h2_dat_compressed$Scaled)
h2_dat_compressed <- h2_dat_compressed[h2_dat_compressed$Ancestry != "EAS",]
FigS2_Part1 <- h2_dat_compressed[h2_dat_compressed$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01"),]

h2_dat <- read.csv("/data/williamsjacr/UKB_WES_Simulation/Sim_Characteristics_NewProp_PropRareVariants.csv")

h2_dat_compressed <- aggregate(.~Ancestry + causal_prop + scaled,data = h2_dat[,c("Ancestry","causal_prop","scaled","h2_rare","Average_Burden_rare")],mean)
colnames(h2_dat_compressed) <- c("Ancestry","Causal_Prop","Scaled","Average_h2","Average_Burden")
h2_dat_compressed_se <- aggregate(h2_rare~Ancestry + causal_prop + scaled,data = h2_dat,function(x){quantile(x,0.025)})
colnames(h2_dat_compressed_se) <- c("Ancestry","Causal_Prop","Scaled","Q_025")
h2_dat_compressed <- inner_join(h2_dat_compressed,h2_dat_compressed_se)
h2_dat_compressed_se <- aggregate(h2_rare~Ancestry + causal_prop + scaled,data = h2_dat,function(x){quantile(x,0.975)})
colnames(h2_dat_compressed_se) <- c("Ancestry","Causal_Prop","Scaled","Q_975")
h2_dat_compressed <- inner_join(h2_dat_compressed,h2_dat_compressed_se)


h2_dat_compressed$Causal_Prop <- paste0("Causal Prop. ",h2_dat_compressed$Causal_Prop)
h2_dat_compressed$Scaled <- paste0("Scaled: ",h2_dat_compressed$Scaled)
h2_dat_compressed <- h2_dat_compressed[h2_dat_compressed$Ancestry != "EAS",]
FigS2_Part2 <- h2_dat_compressed[h2_dat_compressed$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01"),]

#################################################################
### Figure 4 + Supplementary Figure 7c
#################################################################

full_results <- NULL
full_results_Boot <- NULL
full_results_Boot_Comparison <- NULL

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
  full_results_Boot_Comparison <- rbind(full_results_Boot_Comparison,read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/",trait,"_Comparison_Bootstraps.csv")))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results_Boot_Comparison$trait[full_results_Boot_Comparison$trait == "logTG"] <- "log(TG)"
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

lower_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1)] <- paste0(colnames(lower_95)[-c(1)],"_Lower")
upper_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1)] <- paste0(colnames(upper_95)[-c(1)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CT_Lower_95 = c(CI_95$R2_raw_EUR_RICE_vs_CT_Lower,CI_95$R2_raw_SAS_RICE_vs_CT_Lower,CI_95$R2_raw_AFR_RICE_vs_CT_Lower,CI_95$R2_raw_AMR_RICE_vs_CT_Lower),
                    R2_raw_RICE_vs_CT_Upper_95 = c(CI_95$R2_raw_EUR_RICE_vs_CT_Upper,CI_95$R2_raw_SAS_RICE_vs_CT_Upper,CI_95$R2_raw_AFR_RICE_vs_CT_Upper,CI_95$R2_raw_AMR_RICE_vs_CT_Upper),
                    R2_raw_RICE_vs_LDpred2_Lower_95 = c(CI_95$R2_raw_EUR_RICE_vs_LDpred2_Lower,CI_95$R2_raw_SAS_RICE_vs_LDpred2_Lower,CI_95$R2_raw_AFR_RICE_vs_LDpred2_Lower,CI_95$R2_raw_AMR_RICE_vs_LDpred2_Lower),
                    R2_raw_RICE_vs_LDpred2_Upper_95 = c(CI_95$R2_raw_EUR_RICE_vs_LDpred2_Upper,CI_95$R2_raw_SAS_RICE_vs_LDpred2_Upper,CI_95$R2_raw_AFR_RICE_vs_LDpred2_Upper,CI_95$R2_raw_AMR_RICE_vs_LDpred2_Upper),
                    R2_raw_RICE_vs_Lassosum2_Lower_95 = c(CI_95$R2_raw_EUR_RICE_vs_Lassosum2_Lower,CI_95$R2_raw_SAS_RICE_vs_Lassosum2_Lower,CI_95$R2_raw_AFR_RICE_vs_Lassosum2_Lower,CI_95$R2_raw_AMR_RICE_vs_Lassosum2_Lower),
                    R2_raw_RICE_vs_Lassosum2_Upper_95 = c(CI_95$R2_raw_EUR_RICE_vs_Lassosum2_Upper,CI_95$R2_raw_SAS_RICE_vs_Lassosum2_Upper,CI_95$R2_raw_AFR_RICE_vs_Lassosum2_Upper,CI_95$R2_raw_AMR_RICE_vs_Lassosum2_Upper),
                    R2_adjusted_RICE_vs_CT_Lower_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_CT_Lower,CI_95$R2_adjusted_SAS_RICE_vs_CT_Lower,CI_95$R2_adjusted_AFR_RICE_vs_CT_Lower,CI_95$R2_adjusted_AMR_RICE_vs_CT_Lower),
                    R2_adjusted_RICE_vs_CT_Upper_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_CT_Upper,CI_95$R2_adjusted_SAS_RICE_vs_CT_Upper,CI_95$R2_adjusted_AFR_RICE_vs_CT_Upper,CI_95$R2_adjusted_AMR_RICE_vs_CT_Upper),
                    R2_adjusted_RICE_vs_LDpred2_Lower_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_95$R2_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_95$R2_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_95$R2_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    R2_adjusted_RICE_vs_LDpred2_Upper_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_95$R2_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_95$R2_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_95$R2_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    R2_adjusted_RICE_vs_Lassosum2_Lower_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_95$R2_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_95$R2_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_95$R2_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    R2_adjusted_RICE_vs_Lassosum2_Upper_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_95$R2_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_95$R2_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_95$R2_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_raw_Lower_99 = c(CI_99$beta_raw_EUR_boot_Lower,CI_99$beta_raw_SAS_boot_Lower,CI_99$beta_raw_AFR_boot_Lower,CI_99$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_99 = c(CI_99$beta_raw_EUR_boot_Upper,CI_99$beta_raw_SAS_boot_Upper,CI_99$beta_raw_AFR_boot_Upper,CI_99$beta_raw_AMR_boot_Upper),
                    R2_raw_Lower_99 = c(CI_99$R2_raw_EUR_boot_Lower,CI_99$R2_raw_SAS_boot_Lower,CI_99$R2_raw_AFR_boot_Lower,CI_99$R2_raw_AMR_boot_Lower),
                    R2_raw_Upper_99 = c(CI_99$R2_raw_EUR_boot_Upper,CI_99$R2_raw_SAS_boot_Upper,CI_99$R2_raw_AFR_boot_Upper,CI_99$R2_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower,CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper,CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper),
                    R2_adjusted_Lower_99 = c(CI_99$R2_adjusted_EUR_boot_Lower,CI_99$R2_adjusted_SAS_boot_Lower,CI_99$R2_adjusted_AFR_boot_Lower,CI_99$R2_adjusted_AMR_boot_Lower),
                    R2_adjusted_Upper_99 = c(CI_99$R2_adjusted_EUR_boot_Upper,CI_99$R2_adjusted_SAS_boot_Upper,CI_99$R2_adjusted_AFR_boot_Upper,CI_99$R2_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

lower_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_99)[-c(1)] <- paste0(colnames(lower_99)[-c(1)],"_Lower")
upper_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_99)[-c(1)] <- paste0(colnames(upper_99)[-c(1)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CT_Lower_99 = c(CI_99$R2_raw_EUR_RICE_vs_CT_Lower,CI_99$R2_raw_SAS_RICE_vs_CT_Lower,CI_99$R2_raw_AFR_RICE_vs_CT_Lower,CI_99$R2_raw_AMR_RICE_vs_CT_Lower),
                    R2_raw_RICE_vs_CT_Upper_99 = c(CI_99$R2_raw_EUR_RICE_vs_CT_Upper,CI_99$R2_raw_SAS_RICE_vs_CT_Upper,CI_99$R2_raw_AFR_RICE_vs_CT_Upper,CI_99$R2_raw_AMR_RICE_vs_CT_Upper),
                    R2_raw_RICE_vs_LDpred2_Lower_99 = c(CI_99$R2_raw_EUR_RICE_vs_LDpred2_Lower,CI_99$R2_raw_SAS_RICE_vs_LDpred2_Lower,CI_99$R2_raw_AFR_RICE_vs_LDpred2_Lower,CI_99$R2_raw_AMR_RICE_vs_LDpred2_Lower),
                    R2_raw_RICE_vs_LDpred2_Upper_99 = c(CI_99$R2_raw_EUR_RICE_vs_LDpred2_Upper,CI_99$R2_raw_SAS_RICE_vs_LDpred2_Upper,CI_99$R2_raw_AFR_RICE_vs_LDpred2_Upper,CI_99$R2_raw_AMR_RICE_vs_LDpred2_Upper),
                    R2_raw_RICE_vs_Lassosum2_Lower_99 = c(CI_99$R2_raw_EUR_RICE_vs_Lassosum2_Lower,CI_99$R2_raw_SAS_RICE_vs_Lassosum2_Lower,CI_99$R2_raw_AFR_RICE_vs_Lassosum2_Lower,CI_99$R2_raw_AMR_RICE_vs_Lassosum2_Lower),
                    R2_raw_RICE_vs_Lassosum2_Upper_99 = c(CI_99$R2_raw_EUR_RICE_vs_Lassosum2_Upper,CI_99$R2_raw_SAS_RICE_vs_Lassosum2_Upper,CI_99$R2_raw_AFR_RICE_vs_Lassosum2_Upper,CI_99$R2_raw_AMR_RICE_vs_Lassosum2_Upper),
                    R2_adjusted_RICE_vs_CT_Lower_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_CT_Lower,CI_99$R2_adjusted_SAS_RICE_vs_CT_Lower,CI_99$R2_adjusted_AFR_RICE_vs_CT_Lower,CI_99$R2_adjusted_AMR_RICE_vs_CT_Lower),
                    R2_adjusted_RICE_vs_CT_Upper_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_CT_Upper,CI_99$R2_adjusted_SAS_RICE_vs_CT_Upper,CI_99$R2_adjusted_AFR_RICE_vs_CT_Upper,CI_99$R2_adjusted_AMR_RICE_vs_CT_Upper),
                    R2_adjusted_RICE_vs_LDpred2_Lower_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_99$R2_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_99$R2_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_99$R2_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    R2_adjusted_RICE_vs_LDpred2_Upper_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_99$R2_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_99$R2_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_99$R2_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    R2_adjusted_RICE_vs_Lassosum2_Lower_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_99$R2_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_99$R2_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_99$R2_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    R2_adjusted_RICE_vs_Lassosum2_Upper_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_99$R2_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_99$R2_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_99$R2_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_99)

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower_95 = full_results$beta_raw_Lower_95, upper_95 = full_results$beta_raw_Upper_95,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower_95 = full_results$beta_adjusted_Lower_95, upper_95 = full_results$beta_adjusted_Upper_95,method = full_results$Method,Standardization = "Using PCs 1-5"))

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-CV"
full_results$group2 <- "RICE-CV"
full_results$p.signif_beta <- ""
full_results$p.signif_beta[full_results$Method == "RICE-CV"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV"] > 0,"**",""))
full_results$position <- NA
full_results$position[full_results$Method == "RICE-CV"] <- full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"] + 0.03
ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"],full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

Fig4_Betas <- full_results

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV",]
full_results$Method[full_results$Method == "RICE-CV"] <- "RICE"
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

full_results$group1 <- "RICE"
full_results$group2 <- "RICE"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta2 <- ""

for(trait in c("BMI","Height","LDL","log(TG)","TC","HDL")){
  for(anc in c("AFR","EUR","SAS","AMR")){
    tmp <- full_results[full_results$ancestry == anc & full_results$trait == trait,]
    max_R2_notRICE <- max(tmp$R2_adjusted[tmp$Method != "RICE"])
    Best_Method <- tmp$Method[tmp$R2_adjusted == max_R2_notRICE]
    Improved_R2 <- round((tmp$R2_adjusted[tmp$Method == "RICE"]/max_R2_notRICE - 1)*100,digits = 2)
    
    if(Best_Method == "CT"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0("(***)"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0("(**)"),""))
    }else if(Best_Method == "LDpred2"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else{
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }
  }
}

full_results$position1 <- NA
full_results$position2 <- NA
full_results$position1[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.01
full_results$position2[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.07
ylim <- max(c(full_results$R2_adjusted)) + 0.08


full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

FigS7_Continuous_R2 <- full_results

#################################################################
### Supplementary Figures 7a + 7b
#################################################################

full_results <- NULL
full_results_Boot <- NULL
full_results_Boot_Comparison <- NULL

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
  full_results_Boot_Comparison <- rbind(full_results_Boot_Comparison,read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/",trait,"_Comparison_Bootstraps.csv")))
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

lower_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1)] <- paste0(colnames(lower_95)[-c(1)],"_Lower")
upper_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1)] <- paste0(colnames(upper_95)[-c(1)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = "RICE-CV",
                    AUC_raw_RICE_vs_CT_Lower_95 = c(CI_95$AUC_raw_EUR_RICE_vs_CT_Lower,CI_95$AUC_raw_SAS_RICE_vs_CT_Lower,CI_95$AUC_raw_AFR_RICE_vs_CT_Lower,CI_95$AUC_raw_AMR_RICE_vs_CT_Lower),
                    AUC_raw_RICE_vs_CT_Upper_95 = c(CI_95$AUC_raw_EUR_RICE_vs_CT_Upper,CI_95$AUC_raw_SAS_RICE_vs_CT_Upper,CI_95$AUC_raw_AFR_RICE_vs_CT_Upper,CI_95$AUC_raw_AMR_RICE_vs_CT_Upper),
                    AUC_raw_RICE_vs_LDpred2_Lower_95 = c(CI_95$AUC_raw_EUR_RICE_vs_LDpred2_Lower,CI_95$AUC_raw_SAS_RICE_vs_LDpred2_Lower,CI_95$AUC_raw_AFR_RICE_vs_LDpred2_Lower,CI_95$AUC_raw_AMR_RICE_vs_LDpred2_Lower),
                    AUC_raw_RICE_vs_LDpred2_Upper_95 = c(CI_95$AUC_raw_EUR_RICE_vs_LDpred2_Upper,CI_95$AUC_raw_SAS_RICE_vs_LDpred2_Upper,CI_95$AUC_raw_AFR_RICE_vs_LDpred2_Upper,CI_95$AUC_raw_AMR_RICE_vs_LDpred2_Upper),
                    AUC_raw_RICE_vs_Lassosum2_Lower_95 = c(CI_95$AUC_raw_EUR_RICE_vs_Lassosum2_Lower,CI_95$AUC_raw_SAS_RICE_vs_Lassosum2_Lower,CI_95$AUC_raw_AFR_RICE_vs_Lassosum2_Lower,CI_95$AUC_raw_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_raw_RICE_vs_Lassosum2_Upper_95 = c(CI_95$AUC_raw_EUR_RICE_vs_Lassosum2_Upper,CI_95$AUC_raw_SAS_RICE_vs_Lassosum2_Upper,CI_95$AUC_raw_AFR_RICE_vs_Lassosum2_Upper,CI_95$AUC_raw_AMR_RICE_vs_Lassosum2_Upper),
                    AUC_adjusted_RICE_vs_CT_Lower_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_CT_Lower,CI_95$AUC_adjusted_SAS_RICE_vs_CT_Lower,CI_95$AUC_adjusted_AFR_RICE_vs_CT_Lower,CI_95$AUC_adjusted_AMR_RICE_vs_CT_Lower),
                    AUC_adjusted_RICE_vs_CT_Upper_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_CT_Upper,CI_95$AUC_adjusted_SAS_RICE_vs_CT_Upper,CI_95$AUC_adjusted_AFR_RICE_vs_CT_Upper,CI_95$AUC_adjusted_AMR_RICE_vs_CT_Upper),
                    AUC_adjusted_RICE_vs_LDpred2_Lower_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_95$AUC_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_95$AUC_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_95$AUC_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    AUC_adjusted_RICE_vs_LDpred2_Upper_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_95$AUC_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_95$AUC_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_95$AUC_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    AUC_adjusted_RICE_vs_Lassosum2_Lower_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_95$AUC_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_95$AUC_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_95$AUC_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_adjusted_RICE_vs_Lassosum2_Upper_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_95$AUC_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_95$AUC_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_95$AUC_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_raw_Lower_99 = c(CI_99$beta_raw_EUR_boot_Lower,CI_99$beta_raw_SAS_boot_Lower,CI_99$beta_raw_AFR_boot_Lower,CI_99$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_99 = c(CI_99$beta_raw_EUR_boot_Upper,CI_99$beta_raw_SAS_boot_Upper,CI_99$beta_raw_AFR_boot_Upper,CI_99$beta_raw_AMR_boot_Upper),
                    AUC_raw_Lower_99 = c(CI_99$AUC_raw_EUR_boot_Lower,CI_99$AUC_raw_SAS_boot_Lower,CI_99$AUC_raw_AFR_boot_Lower,CI_99$AUC_raw_AMR_boot_Lower),
                    AUC_raw_Upper_99 = c(CI_99$AUC_raw_EUR_boot_Upper,CI_99$AUC_raw_SAS_boot_Upper,CI_99$AUC_raw_AFR_boot_Upper,CI_99$AUC_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower,CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper,CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper),
                    AUC_adjusted_Lower_99 = c(CI_99$AUC_adjusted_EUR_boot_Lower,CI_99$AUC_adjusted_SAS_boot_Lower,CI_99$AUC_adjusted_AFR_boot_Lower,CI_99$AUC_adjusted_AMR_boot_Lower),
                    AUC_adjusted_Upper_99 = c(CI_99$AUC_adjusted_EUR_boot_Upper,CI_99$AUC_adjusted_SAS_boot_Upper,CI_99$AUC_adjusted_AFR_boot_Upper,CI_99$AUC_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

lower_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_99)[-c(1)] <- paste0(colnames(lower_99)[-c(1)],"_Lower")
upper_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_99)[-c(1)] <- paste0(colnames(upper_99)[-c(1)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = "RICE-CV",
                    AUC_raw_RICE_vs_CT_Lower_99 = c(CI_99$AUC_raw_EUR_RICE_vs_CT_Lower,CI_99$AUC_raw_SAS_RICE_vs_CT_Lower,CI_99$AUC_raw_AFR_RICE_vs_CT_Lower,CI_99$AUC_raw_AMR_RICE_vs_CT_Lower),
                    AUC_raw_RICE_vs_CT_Upper_99 = c(CI_99$AUC_raw_EUR_RICE_vs_CT_Upper,CI_99$AUC_raw_SAS_RICE_vs_CT_Upper,CI_99$AUC_raw_AFR_RICE_vs_CT_Upper,CI_99$AUC_raw_AMR_RICE_vs_CT_Upper),
                    AUC_raw_RICE_vs_LDpred2_Lower_99 = c(CI_99$AUC_raw_EUR_RICE_vs_LDpred2_Lower,CI_99$AUC_raw_SAS_RICE_vs_LDpred2_Lower,CI_99$AUC_raw_AFR_RICE_vs_LDpred2_Lower,CI_99$AUC_raw_AMR_RICE_vs_LDpred2_Lower),
                    AUC_raw_RICE_vs_LDpred2_Upper_99 = c(CI_99$AUC_raw_EUR_RICE_vs_LDpred2_Upper,CI_99$AUC_raw_SAS_RICE_vs_LDpred2_Upper,CI_99$AUC_raw_AFR_RICE_vs_LDpred2_Upper,CI_99$AUC_raw_AMR_RICE_vs_LDpred2_Upper),
                    AUC_raw_RICE_vs_Lassosum2_Lower_99 = c(CI_99$AUC_raw_EUR_RICE_vs_Lassosum2_Lower,CI_99$AUC_raw_SAS_RICE_vs_Lassosum2_Lower,CI_99$AUC_raw_AFR_RICE_vs_Lassosum2_Lower,CI_99$AUC_raw_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_raw_RICE_vs_Lassosum2_Upper_99 = c(CI_99$AUC_raw_EUR_RICE_vs_Lassosum2_Upper,CI_99$AUC_raw_SAS_RICE_vs_Lassosum2_Upper,CI_99$AUC_raw_AFR_RICE_vs_Lassosum2_Upper,CI_99$AUC_raw_AMR_RICE_vs_Lassosum2_Upper),
                    AUC_adjusted_RICE_vs_CT_Lower_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_CT_Lower,CI_99$AUC_adjusted_SAS_RICE_vs_CT_Lower,CI_99$AUC_adjusted_AFR_RICE_vs_CT_Lower,CI_99$AUC_adjusted_AMR_RICE_vs_CT_Lower),
                    AUC_adjusted_RICE_vs_CT_Upper_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_CT_Upper,CI_99$AUC_adjusted_SAS_RICE_vs_CT_Upper,CI_99$AUC_adjusted_AFR_RICE_vs_CT_Upper,CI_99$AUC_adjusted_AMR_RICE_vs_CT_Upper),
                    AUC_adjusted_RICE_vs_LDpred2_Lower_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_99$AUC_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_99$AUC_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_99$AUC_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    AUC_adjusted_RICE_vs_LDpred2_Upper_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_99$AUC_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_99$AUC_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_99$AUC_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    AUC_adjusted_RICE_vs_Lassosum2_Lower_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_99$AUC_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_99$AUC_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_99$AUC_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_adjusted_RICE_vs_Lassosum2_Upper_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_99$AUC_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_99$AUC_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_99$AUC_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_99)

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower_95 = full_results$beta_raw_Lower_95, upper_95 = full_results$beta_raw_Upper_95,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower_95 = full_results$beta_adjusted_Lower_95, upper_95 = full_results$beta_adjusted_Upper_95,method = full_results$Method,Standardization = "Using PCs 1-5"))


full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-CV"
full_results$group2 <- "RICE-CV"
full_results$p.signif_beta <- ""
full_results$p.signif_beta[full_results$Method == "RICE-CV"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV"] > 0,"**",""))
full_results$position <- NA
full_results$position[full_results$Method == "RICE-CV"] <- full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"] + 0.03
ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"],full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

FigS7_Binary_Betas <- full_results

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV",]
full_results$Method[full_results$Method == "RICE-CV"] <- "RICE"
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

full_results$group1 <- "RICE"
full_results$group2 <- "RICE"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta2 <- ""

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  for(anc in c("AFR","EUR","SAS","AMR")){
    tmp <- full_results[full_results$ancestry == anc & full_results$trait == trait,]
    max_AUC_notRICE <- max(tmp$AUC_adjusted[tmp$Method != "RICE"])
    Best_Method <- tmp$Method[tmp$AUC_adjusted == max_AUC_notRICE]
    Improved_AUC <- round((tmp$AUC_adjusted[tmp$Method == "RICE"]/max_AUC_notRICE - 1)*100,digits = 2)
    
    if(Best_Method == "CT"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0("(***)"),ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0("(**)"),""))
    }else if(Best_Method == "LDpred2"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else{
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }
  }
}

full_results$position1 <- NA
full_results$position2 <- NA
full_results$position1[full_results$Method == "RICE"] <- full_results$AUC_adjusted[full_results$Method == "RICE"] + 0.01
full_results$position2[full_results$Method == "RICE"] <- full_results$AUC_adjusted[full_results$Method == "RICE"] + 0.07
ylim <- max(c(full_results$AUC_adjusted)) + 0.08

full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

FigS7_Binary_AUC <- full_results

#################################################################
### Figure 5 + Supplementary Figure 9
#################################################################

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
Fig5_PlotA <- NULL
Fig5_PlotB <- NULL
Fig5_PlotC <- NULL
Fig5_StatC <- NULL
for(trait in continuous_traits){
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
  pheno_validation$y_validation <- NA
  pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
  
  pheno_validation <- pheno_validation[!is.na(pheno_validation[,trait]),]
  
  CV_RV_PRS_raw <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for(i in c("RV_PRS","CV_PRS")){
    tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      CV_RV_PRS_adjusted[,i] <- 0
    }else{
      CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
    }
  }
  
  NRI_Data_Continuous <- NULL
  
  for(risk in c(0.05,0.1)){
    
    truth_HighRisk <- which(pheno_validation$y_validation > quantile(pheno_validation$y_validation,0.9))
    CV_HighRisk_RV_HighRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_HighRisk_RV_NotHighRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,1 - risk)))
    
    truth_NotHighRisk <- which(pheno_validation$y_validation < quantile(pheno_validation$y_validation,0.9))
    CV_NotHighRisk_RV_HighRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_NotHighRisk_RV_NotHighRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,1 - risk)))
    
    tmp <- data.frame(trait = trait, risk = paste0(100*risk,"%"),
                      
                      Percent_A = 100*length(CV_NotHighRisk_RV_NotHighRisk)/nrow(pheno_validation),
                      Mean_A = mean(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk]),
                      SE_A = sd(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk])/sqrt(length(CV_NotHighRisk_RV_NotHighRisk)),
                      Total_Capture_A = sum(CV_NotHighRisk_RV_NotHighRisk %in% truth_HighRisk),
                      
                      Percent_B = 100*length(CV_HighRisk_RV_NotHighRisk)/nrow(pheno_validation),
                      Mean_B = mean(pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk]),
                      SE_B = sd(pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk])/sqrt(length(CV_HighRisk_RV_NotHighRisk)),
                      Total_Capture_B = sum(CV_HighRisk_RV_NotHighRisk %in% truth_HighRisk),
                      
                      Percent_C = 100*length(CV_NotHighRisk_RV_HighRisk)/nrow(pheno_validation),
                      Mean_C = mean(pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk]),
                      SE_C = sd(pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk])/sqrt(length(CV_NotHighRisk_RV_HighRisk)),
                      Total_Capture_C = sum(CV_NotHighRisk_RV_HighRisk %in% truth_HighRisk),
                      
                      Percent_D = 100*length(CV_HighRisk_RV_HighRisk)/nrow(pheno_validation),
                      Mean_D = mean(pheno_validation$y_validation[CV_HighRisk_RV_HighRisk]),
                      SE_D = sd(pheno_validation$y_validation[CV_HighRisk_RV_HighRisk])/sqrt(length(CV_HighRisk_RV_HighRisk)),
                      Total_Capture_D = sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk),
                      
                      Mean_BD = mean(pheno_validation$y_validation[c(CV_HighRisk_RV_HighRisk,CV_HighRisk_RV_NotHighRisk)]))
    
    tmp$C_minus_A <- (tmp$Mean_C - tmp$Mean_A)/sd(pheno_validation$y_validation)
    tmp$B_minus_A <- (tmp$Mean_B - tmp$Mean_A)/sd(pheno_validation$y_validation)
    tmp$BD_minus_A <- (tmp$Mean_BD - tmp$Mean_A)/sd(pheno_validation$y_validation)
    
    tmp$A_vs_B <- t.test(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk], pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk], alternative = "two.sided", var.equal = FALSE)$p.value
    tmp$A_vs_C <- t.test(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk], pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk], alternative = "two.sided", var.equal = FALSE)$p.value
    tmp$A_vs_D <- t.test(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk], pheno_validation$y_validation[CV_HighRisk_RV_HighRisk], alternative = "two.sided", var.equal = FALSE)$p.value
    
    NRI_Data_Continuous <- rbind(NRI_Data_Continuous,tmp)
  } 
  
  CV_RV_PRS_adjusted$Common_Bin <- 9
  CV_RV_PRS_adjusted$Rare_Bin <- 9
  
  Common_quants <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(0,.1,.2,.3,.4,.6,.7,.8,.9))
  
  for(i in 1:8){
    CV_RV_PRS_adjusted$Common_Bin[CV_RV_PRS_adjusted$CV_PRS < unname(Common_quants[i + 1]) & CV_RV_PRS_adjusted$CV_PRS >= unname(Common_quants[i])] <- i
  }
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted$y_validation <- scale(CV_RV_PRS_adjusted$y_validation)
  
  CV_RV_PRS_adjusted <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted$Rare_Bin <- factor(CV_RV_PRS_adjusted$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  CV_RV_PRS_adjusted_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted,mean)
  
  colnames(CV_RV_PRS_adjusted) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted <- inner_join(CV_RV_PRS_adjusted,CV_RV_PRS_adjusted_se)
  
  colnames(CV_RV_PRS_adjusted) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  
  ymin <- round(min(c(CV_RV_PRS_adjusted$Mean - CV_RV_PRS_adjusted$SE)) - 0.05,2)
  ymax <- round(max(c(CV_RV_PRS_adjusted$Mean + CV_RV_PRS_adjusted$SE)) + 0.05,2)
  
  Fig5_PlotA <- rbind(Fig5_PlotA,cbind(CV_RV_PRS_adjusted,trait))
  
  plot_data <- data.frame(Method = rep(c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
                          trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
                          risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
                          value = c(NRI_Data_Continuous$Total_Capture_A,NRI_Data_Continuous$Total_Capture_B,NRI_Data_Continuous$Total_Capture_C,NRI_Data_Continuous$Total_Capture_D))
  
  plot_data$Method <- factor(plot_data$Method, levels = c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"))
  
  plot_data_sub <- plot_data[plot_data$trait == trait & plot_data$risk == "5%", ] %>%
    arrange(Method) %>%
    mutate(
      fraction = value / sum(value),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, -1)),
      label_pos = (ymin + ymax) / 2,
      label = paste0(round(100*value/sum(value), 1), "%")
    )
  
  plot_data_sub$Method <- factor(plot_data_sub$Method, levels = c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"))
  
  Fig5_PlotB <- rbind(Fig5_PlotB,cbind(plot_data_sub,trait))
  
  plot_data <- data.frame(Method = rep(c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
                          trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
                          risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
                          mean = c(NRI_Data_Continuous$Mean_A,NRI_Data_Continuous$Mean_B,NRI_Data_Continuous$Mean_C,NRI_Data_Continuous$Mean_D),
                          se = c(NRI_Data_Continuous$SE_A,NRI_Data_Continuous$SE_B,NRI_Data_Continuous$SE_C,NRI_Data_Continuous$SE_D))
  
  plot_data$Method <- factor(plot_data$Method, levels = c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"))
  
  if(trait %in% c("LDL","BMI","TC")){
    stat.test <- data.frame(
      group1 = "Low CV PRS, Low RV PRS",
      group2 = rep(c("High CV PRS, Low RV PRS", "Low CV PRS, High RV PRS", "High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
      trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
      risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
      p.adj = signif(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D),3),      # Adjusted p-values
      y.position = as.vector(outer(apply(cbind(NRI_Data_Continuous$Mean_B,NRI_Data_Continuous$Mean_C,NRI_Data_Continuous$Mean_D),1,max),c(1.3,1.4,1.5),"*")),     # Height of brackets
      p.adj.signif = ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) > 0.05, "",ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) < 0.01,"**","*"))  # Significance symbols
    ) 
  }else{
    stat.test <- data.frame(
      group1 = "Low CV PRS, Low RV PRS",
      group2 = rep(c("High CV PRS, Low RV PRS", "Low CV PRS, High RV PRS", "High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
      trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
      risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
      p.adj = signif(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D),3),      # Adjusted p-values
      y.position = as.vector(outer(apply(cbind(NRI_Data_Continuous$Mean_B,NRI_Data_Continuous$Mean_C,NRI_Data_Continuous$Mean_D),1,max),c(1.1,1.2,1.3),"*")),     # Height of brackets
      p.adj.signif = ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) > 0.05, "",ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) < 0.01,"**","*"))  # Significance symbols
    )
  }
  
  Fig5_PlotC <- rbind(Fig5_PlotC,plot_data[plot_data$trait == trait & plot_data$risk == "5%",])
  Fig5_StatC <- rbind(Fig5_StatC,stat.test[stat.test$trait == trait & stat.test$risk == "5%",])
}

FigS9_PlotA <- Fig5_PlotA[Fig5_PlotA$trait != "HDL",]
Fig5_PlotA <- Fig5_PlotA[Fig5_PlotA$trait == "HDL",]

FigS9_PlotB <- Fig5_PlotB[Fig5_PlotB$trait != "HDL",]
Fig5_PlotB <- Fig5_PlotB[Fig5_PlotB$trait == "HDL",]

FigS9_PlotC <- Fig5_PlotC[Fig5_PlotC$trait != "HDL",]
Fig5_PlotC <- Fig5_PlotC[Fig5_PlotC$trait == "HDL",]

FigS9_StatC <- Fig5_StatC[Fig5_StatC$trait != "HDL",]
Fig5_StatC <- Fig5_StatC[Fig5_StatC$trait == "HDL",]


#################################################################
### Supplementary Figure 8
#################################################################

FigS8_EUR <- NULL
FigS8_AMR <- NULL
FigS8_AFR <- NULL
FigS8_SAS <- NULL

for(trait in c("BMI","HDL","LDL","logTG","TC","Height")){
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
  pheno_validation$y_validation <- NA
  pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
  
  CV_RV_PRS_raw <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for(i in c("RV_PRS","CV_PRS")){
    tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      CV_RV_PRS_adjusted[,i] <- 0
    }else{
      CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
    }
  }
  
  
  CV_RV_PRS_adjusted$Common_Bin <- 9
  CV_RV_PRS_adjusted$Rare_Bin <- 9
  
  Common_quants <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(0,.1,.2,.3,.4,.6,.7,.8,.9))
  
  for(i in 1:8){
    CV_RV_PRS_adjusted$Common_Bin[CV_RV_PRS_adjusted$CV_PRS < unname(Common_quants[i + 1]) & CV_RV_PRS_adjusted$CV_PRS >= unname(Common_quants[i])] <- i
  }
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_EUR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_EUR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_EUR$y_validation <- scale(CV_RV_PRS_adjusted_EUR$y_validation)
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted_EUR[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin <- factor(CV_RV_PRS_adjusted_EUR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_AMR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AMR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_AMR$y_validation <- scale(CV_RV_PRS_adjusted_AMR$y_validation)
  
  CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted_AMR[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_AMR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AMR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_AFR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AFR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_AFR$y_validation <- scale(CV_RV_PRS_adjusted_AFR$y_validation)
  
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted_AFR[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_AFR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AFR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_SAS$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_SAS$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_SAS$y_validation <- scale(CV_RV_PRS_adjusted_SAS$y_validation)
  
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted_SAS[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_SAS$Rare_Bin <- factor(CV_RV_PRS_adjusted_SAS$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  CV_RV_PRS_adjusted_EUR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_EUR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_EUR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_EUR,mean)
  
  colnames(CV_RV_PRS_adjusted_EUR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_EUR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_EUR <- inner_join(CV_RV_PRS_adjusted_EUR,CV_RV_PRS_adjusted_EUR_se)
  
  CV_RV_PRS_adjusted_AMR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AMR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_AMR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AMR,mean)
  
  colnames(CV_RV_PRS_adjusted_AMR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_AMR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_AMR <- inner_join(CV_RV_PRS_adjusted_AMR,CV_RV_PRS_adjusted_AMR_se)
  
  CV_RV_PRS_adjusted_AFR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AFR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_AFR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AFR,mean)
  
  colnames(CV_RV_PRS_adjusted_AFR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_AFR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_AFR <- inner_join(CV_RV_PRS_adjusted_AFR,CV_RV_PRS_adjusted_AFR_se)
  
  CV_RV_PRS_adjusted_SAS_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_SAS,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_SAS <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_SAS,mean)
  
  colnames(CV_RV_PRS_adjusted_SAS) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_SAS_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_SAS <- inner_join(CV_RV_PRS_adjusted_SAS,CV_RV_PRS_adjusted_SAS_se)
  
  
  colnames(CV_RV_PRS_adjusted_EUR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_AMR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_AFR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_SAS) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  
  ymin_EUR <- round(min(c(CV_RV_PRS_adjusted_EUR$Mean - CV_RV_PRS_adjusted_EUR$SE)) - 0.05,2)
  ymax_EUR <- round(max(c(CV_RV_PRS_adjusted_EUR$Mean + CV_RV_PRS_adjusted_EUR$SE)) + 0.05,2)
  
  ymin_AMR <- round(min(c(CV_RV_PRS_adjusted_AMR$Mean - CV_RV_PRS_adjusted_AMR$SE)) - 0.05,2)
  ymax_AMR <- round(max(c(CV_RV_PRS_adjusted_AMR$Mean + CV_RV_PRS_adjusted_AMR$SE)) + 0.05,2)
  
  ymin_AFR <- round(min(c(CV_RV_PRS_adjusted_AFR$Mean - CV_RV_PRS_adjusted_AFR$SE)) - 0.05,2)
  ymax_AFR <- round(max(c(CV_RV_PRS_adjusted_AFR$Mean + CV_RV_PRS_adjusted_AFR$SE)) + 0.05,2)
  
  ymin_SAS <- round(min(c(CV_RV_PRS_adjusted_SAS$Mean - CV_RV_PRS_adjusted_SAS$SE)) - 0.05,2)
  ymax_SAS <- round(max(c(CV_RV_PRS_adjusted_SAS$Mean + CV_RV_PRS_adjusted_SAS$SE)) + 0.05,2)
  
  FigS8_EUR <- rbind(FigS8_EUR,cbind(CV_RV_PRS_adjusted_EUR,trait))
  FigS8_AMR <- rbind(FigS8_AMR,cbind(CV_RV_PRS_adjusted_AMR,trait))
  FigS8_AFR <- rbind(FigS8_AFR,cbind(CV_RV_PRS_adjusted_AFR,trait))
  FigS8_SAS <- rbind(FigS8_SAS,cbind(CV_RV_PRS_adjusted_SAS,trait))
  
}


#################################################################
### Supplementary Figure 10
#################################################################

full_results <- NULL

for(trait in c("TC","HDL","LDL","logTG")){
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_Lipids_Gene_Best_Betas.csv"))
  tmp$Method <- "High-Penetrance Genes"
  tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  full_results <- rbind(full_results,tmp) 
  
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"Best_Betas.csv"))
  tmp$Method <- "RICE-RV"
  tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  full_results <- rbind(full_results,tmp) 
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"

FigS10 <- full_results

#################################################################
### Supplementary Figure 11
#################################################################

full_results <- NULL

options(scipen = 0)

for(thresholds in c(1e-5,1e-4,1e-3,1e-2)){
  for(trait in c("BMI","LDL","HDL","logTG","TC","Height","Breast","Prostate","CAD","T2D","Asthma")){
    tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Threshold_Sensitivity_RICE_RV/",trait,"_",thresholds,"Best_Betas.csv"))
    tmp$Threshold <- thresholds
    tmp <- tmp[,c("Threshold","trait","ancestry","beta_adjusted")]
    full_results <- rbind(full_results,tmp) 
  }
}

options("scipen"=100, "digits"=4)
full_results$Threshold <- as.character(full_results$Threshold)
options(scipen = 0)

full_results$Threshold <- factor(full_results$Threshold,levels = c("0.00001","0.0001","0.001","0.01"))

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"

full_results_continuous <- full_results[full_results$trait %in% c("BMI","Height","HDL","LDL","log(TG)","TC"),]
full_results_binary <- full_results[full_results$trait %in% c("Asthma","CAD","T2D","Breast","Prostate"),]

full_results_continuous$trait <- factor(full_results_continuous$trait,levels = c("BMI","HDL","Height","LDL","log(TG)","TC"))
full_results_binary$trait <- factor(full_results_binary$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))

full_results_continuous$beta_adjusted[full_results_continuous$beta_adjusted < 0] <- 0
full_results_binary$beta_adjusted[full_results_binary$beta_adjusted < 0] <- 0

FigS11 <- full_results_continuous

#################################################################
### Supplementary Figure 12 + 14
#################################################################

full_results <- NULL
full_results_Boot <- NULL
full_results_Boot_Comparison <- NULL

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
  full_results_Boot_Comparison <- rbind(full_results_Boot_Comparison,read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/",trait,"_Comparison_Bootstraps.csv")))
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

lower_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1)] <- paste0(colnames(lower_95)[-c(1)],"_Lower")
upper_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1)] <- paste0(colnames(upper_95)[-c(1)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = "RICE-CV",
                    AUC_raw_RICE_vs_CT_Lower_95 = c(CI_95$AUC_raw_EUR_RICE_vs_CT_Lower,CI_95$AUC_raw_SAS_RICE_vs_CT_Lower,CI_95$AUC_raw_AFR_RICE_vs_CT_Lower,CI_95$AUC_raw_AMR_RICE_vs_CT_Lower),
                    AUC_raw_RICE_vs_CT_Upper_95 = c(CI_95$AUC_raw_EUR_RICE_vs_CT_Upper,CI_95$AUC_raw_SAS_RICE_vs_CT_Upper,CI_95$AUC_raw_AFR_RICE_vs_CT_Upper,CI_95$AUC_raw_AMR_RICE_vs_CT_Upper),
                    AUC_raw_RICE_vs_LDpred2_Lower_95 = c(CI_95$AUC_raw_EUR_RICE_vs_LDpred2_Lower,CI_95$AUC_raw_SAS_RICE_vs_LDpred2_Lower,CI_95$AUC_raw_AFR_RICE_vs_LDpred2_Lower,CI_95$AUC_raw_AMR_RICE_vs_LDpred2_Lower),
                    AUC_raw_RICE_vs_LDpred2_Upper_95 = c(CI_95$AUC_raw_EUR_RICE_vs_LDpred2_Upper,CI_95$AUC_raw_SAS_RICE_vs_LDpred2_Upper,CI_95$AUC_raw_AFR_RICE_vs_LDpred2_Upper,CI_95$AUC_raw_AMR_RICE_vs_LDpred2_Upper),
                    AUC_raw_RICE_vs_Lassosum2_Lower_95 = c(CI_95$AUC_raw_EUR_RICE_vs_Lassosum2_Lower,CI_95$AUC_raw_SAS_RICE_vs_Lassosum2_Lower,CI_95$AUC_raw_AFR_RICE_vs_Lassosum2_Lower,CI_95$AUC_raw_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_raw_RICE_vs_Lassosum2_Upper_95 = c(CI_95$AUC_raw_EUR_RICE_vs_Lassosum2_Upper,CI_95$AUC_raw_SAS_RICE_vs_Lassosum2_Upper,CI_95$AUC_raw_AFR_RICE_vs_Lassosum2_Upper,CI_95$AUC_raw_AMR_RICE_vs_Lassosum2_Upper),
                    AUC_adjusted_RICE_vs_CT_Lower_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_CT_Lower,CI_95$AUC_adjusted_SAS_RICE_vs_CT_Lower,CI_95$AUC_adjusted_AFR_RICE_vs_CT_Lower,CI_95$AUC_adjusted_AMR_RICE_vs_CT_Lower),
                    AUC_adjusted_RICE_vs_CT_Upper_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_CT_Upper,CI_95$AUC_adjusted_SAS_RICE_vs_CT_Upper,CI_95$AUC_adjusted_AFR_RICE_vs_CT_Upper,CI_95$AUC_adjusted_AMR_RICE_vs_CT_Upper),
                    AUC_adjusted_RICE_vs_LDpred2_Lower_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_95$AUC_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_95$AUC_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_95$AUC_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    AUC_adjusted_RICE_vs_LDpred2_Upper_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_95$AUC_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_95$AUC_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_95$AUC_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    AUC_adjusted_RICE_vs_Lassosum2_Lower_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_95$AUC_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_95$AUC_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_95$AUC_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_adjusted_RICE_vs_Lassosum2_Upper_95 = c(CI_95$AUC_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_95$AUC_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_95$AUC_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_95$AUC_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_raw_Lower_99 = c(CI_99$beta_raw_EUR_boot_Lower,CI_99$beta_raw_SAS_boot_Lower,CI_99$beta_raw_AFR_boot_Lower,CI_99$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_99 = c(CI_99$beta_raw_EUR_boot_Upper,CI_99$beta_raw_SAS_boot_Upper,CI_99$beta_raw_AFR_boot_Upper,CI_99$beta_raw_AMR_boot_Upper),
                    AUC_raw_Lower_99 = c(CI_99$AUC_raw_EUR_boot_Lower,CI_99$AUC_raw_SAS_boot_Lower,CI_99$AUC_raw_AFR_boot_Lower,CI_99$AUC_raw_AMR_boot_Lower),
                    AUC_raw_Upper_99 = c(CI_99$AUC_raw_EUR_boot_Upper,CI_99$AUC_raw_SAS_boot_Upper,CI_99$AUC_raw_AFR_boot_Upper,CI_99$AUC_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower,CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper,CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper),
                    AUC_adjusted_Lower_99 = c(CI_99$AUC_adjusted_EUR_boot_Lower,CI_99$AUC_adjusted_SAS_boot_Lower,CI_99$AUC_adjusted_AFR_boot_Lower,CI_99$AUC_adjusted_AMR_boot_Lower),
                    AUC_adjusted_Upper_99 = c(CI_99$AUC_adjusted_EUR_boot_Upper,CI_99$AUC_adjusted_SAS_boot_Upper,CI_99$AUC_adjusted_AFR_boot_Upper,CI_99$AUC_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

lower_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_99)[-c(1)] <- paste0(colnames(lower_99)[-c(1)],"_Lower")
upper_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_99)[-c(1)] <- paste0(colnames(upper_99)[-c(1)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = "RICE-CV",
                    AUC_raw_RICE_vs_CT_Lower_99 = c(CI_99$AUC_raw_EUR_RICE_vs_CT_Lower,CI_99$AUC_raw_SAS_RICE_vs_CT_Lower,CI_99$AUC_raw_AFR_RICE_vs_CT_Lower,CI_99$AUC_raw_AMR_RICE_vs_CT_Lower),
                    AUC_raw_RICE_vs_CT_Upper_99 = c(CI_99$AUC_raw_EUR_RICE_vs_CT_Upper,CI_99$AUC_raw_SAS_RICE_vs_CT_Upper,CI_99$AUC_raw_AFR_RICE_vs_CT_Upper,CI_99$AUC_raw_AMR_RICE_vs_CT_Upper),
                    AUC_raw_RICE_vs_LDpred2_Lower_99 = c(CI_99$AUC_raw_EUR_RICE_vs_LDpred2_Lower,CI_99$AUC_raw_SAS_RICE_vs_LDpred2_Lower,CI_99$AUC_raw_AFR_RICE_vs_LDpred2_Lower,CI_99$AUC_raw_AMR_RICE_vs_LDpred2_Lower),
                    AUC_raw_RICE_vs_LDpred2_Upper_99 = c(CI_99$AUC_raw_EUR_RICE_vs_LDpred2_Upper,CI_99$AUC_raw_SAS_RICE_vs_LDpred2_Upper,CI_99$AUC_raw_AFR_RICE_vs_LDpred2_Upper,CI_99$AUC_raw_AMR_RICE_vs_LDpred2_Upper),
                    AUC_raw_RICE_vs_Lassosum2_Lower_99 = c(CI_99$AUC_raw_EUR_RICE_vs_Lassosum2_Lower,CI_99$AUC_raw_SAS_RICE_vs_Lassosum2_Lower,CI_99$AUC_raw_AFR_RICE_vs_Lassosum2_Lower,CI_99$AUC_raw_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_raw_RICE_vs_Lassosum2_Upper_99 = c(CI_99$AUC_raw_EUR_RICE_vs_Lassosum2_Upper,CI_99$AUC_raw_SAS_RICE_vs_Lassosum2_Upper,CI_99$AUC_raw_AFR_RICE_vs_Lassosum2_Upper,CI_99$AUC_raw_AMR_RICE_vs_Lassosum2_Upper),
                    AUC_adjusted_RICE_vs_CT_Lower_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_CT_Lower,CI_99$AUC_adjusted_SAS_RICE_vs_CT_Lower,CI_99$AUC_adjusted_AFR_RICE_vs_CT_Lower,CI_99$AUC_adjusted_AMR_RICE_vs_CT_Lower),
                    AUC_adjusted_RICE_vs_CT_Upper_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_CT_Upper,CI_99$AUC_adjusted_SAS_RICE_vs_CT_Upper,CI_99$AUC_adjusted_AFR_RICE_vs_CT_Upper,CI_99$AUC_adjusted_AMR_RICE_vs_CT_Upper),
                    AUC_adjusted_RICE_vs_LDpred2_Lower_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_99$AUC_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_99$AUC_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_99$AUC_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    AUC_adjusted_RICE_vs_LDpred2_Upper_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_99$AUC_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_99$AUC_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_99$AUC_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    AUC_adjusted_RICE_vs_Lassosum2_Lower_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_99$AUC_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_99$AUC_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_99$AUC_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    AUC_adjusted_RICE_vs_Lassosum2_Upper_99 = c(CI_99$AUC_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_99$AUC_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_99$AUC_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_99$AUC_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_99)

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower_95 = full_results$beta_raw_Lower_95, upper_95 = full_results$beta_raw_Upper_95,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower_95 = full_results$beta_adjusted_Lower_95, upper_95 = full_results$beta_adjusted_Upper_95,method = full_results$Method,Standardization = "Using PCs 1-5"))

FigS14_Binary <- full_results_stacked


full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-CV"
full_results$group2 <- "RICE-CV"
full_results$p.signif_beta <- ""
full_results$p.signif_beta[full_results$Method == "RICE-CV"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV"] > 0,"**",""))
full_results$position <- NA
full_results$position[full_results$Method == "RICE-CV"] <- full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"] + 0.03
ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"],full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

FigS12_Binary_Betas <- full_results

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV",]
full_results$Method[full_results$Method == "RICE-CV"] <- "RICE"
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

full_results$group1 <- "RICE"
full_results$group2 <- "RICE"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta2 <- ""

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  for(anc in c("AFR","EUR","SAS","AMR")){
    tmp <- full_results[full_results$ancestry == anc & full_results$trait == trait,]
    max_AUC_notRICE <- max(tmp$AUC_adjusted[tmp$Method != "RICE"])
    Best_Method <- tmp$Method[tmp$AUC_adjusted == max_AUC_notRICE]
    Improved_AUC <- round((tmp$AUC_adjusted[tmp$Method == "RICE"]/max_AUC_notRICE - 1)*100,digits = 2)
    
    if(Best_Method == "CT"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0("(***)"),ifelse(tmp$AUC_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0("(**)"),""))
    }else if(Best_Method == "LDpred2"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$AUC_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else{
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_AUC,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$AUC_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }
  }
}

full_results$position1 <- NA
full_results$position2 <- NA
full_results$position1[full_results$Method == "RICE"] <- full_results$AUC_adjusted[full_results$Method == "RICE"] + 0.01
full_results$position2[full_results$Method == "RICE"] <- full_results$AUC_adjusted[full_results$Method == "RICE"] + 0.07
ylim <- max(c(full_results$AUC_adjusted)) + 0.08


full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

FigS12_Binary_AUC <- full_results


full_results <- NULL
full_results_Boot <- NULL
full_results_Boot_Comparison <- NULL

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
  full_results_Boot_Comparison <- rbind(full_results_Boot_Comparison,read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/",trait,"_Comparison_Bootstraps.csv")))
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results_Boot_Comparison$trait[full_results_Boot_Comparison$trait == "logTG"] <- "log(TG)"
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

lower_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1)] <- paste0(colnames(lower_95)[-c(1)],"_Lower")
upper_95 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1)] <- paste0(colnames(upper_95)[-c(1)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_95)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CT_Lower_95 = c(CI_95$R2_raw_EUR_RICE_vs_CT_Lower,CI_95$R2_raw_SAS_RICE_vs_CT_Lower,CI_95$R2_raw_AFR_RICE_vs_CT_Lower,CI_95$R2_raw_AMR_RICE_vs_CT_Lower),
                    R2_raw_RICE_vs_CT_Upper_95 = c(CI_95$R2_raw_EUR_RICE_vs_CT_Upper,CI_95$R2_raw_SAS_RICE_vs_CT_Upper,CI_95$R2_raw_AFR_RICE_vs_CT_Upper,CI_95$R2_raw_AMR_RICE_vs_CT_Upper),
                    R2_raw_RICE_vs_LDpred2_Lower_95 = c(CI_95$R2_raw_EUR_RICE_vs_LDpred2_Lower,CI_95$R2_raw_SAS_RICE_vs_LDpred2_Lower,CI_95$R2_raw_AFR_RICE_vs_LDpred2_Lower,CI_95$R2_raw_AMR_RICE_vs_LDpred2_Lower),
                    R2_raw_RICE_vs_LDpred2_Upper_95 = c(CI_95$R2_raw_EUR_RICE_vs_LDpred2_Upper,CI_95$R2_raw_SAS_RICE_vs_LDpred2_Upper,CI_95$R2_raw_AFR_RICE_vs_LDpred2_Upper,CI_95$R2_raw_AMR_RICE_vs_LDpred2_Upper),
                    R2_raw_RICE_vs_Lassosum2_Lower_95 = c(CI_95$R2_raw_EUR_RICE_vs_Lassosum2_Lower,CI_95$R2_raw_SAS_RICE_vs_Lassosum2_Lower,CI_95$R2_raw_AFR_RICE_vs_Lassosum2_Lower,CI_95$R2_raw_AMR_RICE_vs_Lassosum2_Lower),
                    R2_raw_RICE_vs_Lassosum2_Upper_95 = c(CI_95$R2_raw_EUR_RICE_vs_Lassosum2_Upper,CI_95$R2_raw_SAS_RICE_vs_Lassosum2_Upper,CI_95$R2_raw_AFR_RICE_vs_Lassosum2_Upper,CI_95$R2_raw_AMR_RICE_vs_Lassosum2_Upper),
                    R2_adjusted_RICE_vs_CT_Lower_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_CT_Lower,CI_95$R2_adjusted_SAS_RICE_vs_CT_Lower,CI_95$R2_adjusted_AFR_RICE_vs_CT_Lower,CI_95$R2_adjusted_AMR_RICE_vs_CT_Lower),
                    R2_adjusted_RICE_vs_CT_Upper_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_CT_Upper,CI_95$R2_adjusted_SAS_RICE_vs_CT_Upper,CI_95$R2_adjusted_AFR_RICE_vs_CT_Upper,CI_95$R2_adjusted_AMR_RICE_vs_CT_Upper),
                    R2_adjusted_RICE_vs_LDpred2_Lower_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_95$R2_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_95$R2_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_95$R2_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    R2_adjusted_RICE_vs_LDpred2_Upper_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_95$R2_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_95$R2_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_95$R2_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    R2_adjusted_RICE_vs_Lassosum2_Lower_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_95$R2_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_95$R2_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_95$R2_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    R2_adjusted_RICE_vs_Lassosum2_Upper_95 = c(CI_95$R2_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_95$R2_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_95$R2_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_95$R2_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_raw_Lower_99 = c(CI_99$beta_raw_EUR_boot_Lower,CI_99$beta_raw_SAS_boot_Lower,CI_99$beta_raw_AFR_boot_Lower,CI_99$beta_raw_AMR_boot_Lower),
                    beta_raw_Upper_99 = c(CI_99$beta_raw_EUR_boot_Upper,CI_99$beta_raw_SAS_boot_Upper,CI_99$beta_raw_AFR_boot_Upper,CI_99$beta_raw_AMR_boot_Upper),
                    R2_raw_Lower_99 = c(CI_99$R2_raw_EUR_boot_Lower,CI_99$R2_raw_SAS_boot_Lower,CI_99$R2_raw_AFR_boot_Lower,CI_99$R2_raw_AMR_boot_Lower),
                    R2_raw_Upper_99 = c(CI_99$R2_raw_EUR_boot_Upper,CI_99$R2_raw_SAS_boot_Upper,CI_99$R2_raw_AFR_boot_Upper,CI_99$R2_raw_AMR_boot_Upper),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower,CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper,CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper),
                    R2_adjusted_Lower_99 = c(CI_99$R2_adjusted_EUR_boot_Lower,CI_99$R2_adjusted_SAS_boot_Lower,CI_99$R2_adjusted_AFR_boot_Lower,CI_99$R2_adjusted_AMR_boot_Lower),
                    R2_adjusted_Upper_99 = c(CI_99$R2_adjusted_EUR_boot_Upper,CI_99$R2_adjusted_SAS_boot_Upper,CI_99$R2_adjusted_AFR_boot_Upper,CI_99$R2_adjusted_AMR_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

lower_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.025)})
colnames(lower_99)[-c(1)] <- paste0(colnames(lower_99)[-c(1)],"_Lower")
upper_99 <- aggregate(.~trait,data = full_results_Boot_Comparison,function(x){quantile(x,0.975)})
colnames(upper_99)[-c(1)] <- paste0(colnames(upper_99)[-c(1)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("EUR","SAS","AFR","AMR"),each = nrow(CI_99)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CT_Lower_99 = c(CI_99$R2_raw_EUR_RICE_vs_CT_Lower,CI_99$R2_raw_SAS_RICE_vs_CT_Lower,CI_99$R2_raw_AFR_RICE_vs_CT_Lower,CI_99$R2_raw_AMR_RICE_vs_CT_Lower),
                    R2_raw_RICE_vs_CT_Upper_99 = c(CI_99$R2_raw_EUR_RICE_vs_CT_Upper,CI_99$R2_raw_SAS_RICE_vs_CT_Upper,CI_99$R2_raw_AFR_RICE_vs_CT_Upper,CI_99$R2_raw_AMR_RICE_vs_CT_Upper),
                    R2_raw_RICE_vs_LDpred2_Lower_99 = c(CI_99$R2_raw_EUR_RICE_vs_LDpred2_Lower,CI_99$R2_raw_SAS_RICE_vs_LDpred2_Lower,CI_99$R2_raw_AFR_RICE_vs_LDpred2_Lower,CI_99$R2_raw_AMR_RICE_vs_LDpred2_Lower),
                    R2_raw_RICE_vs_LDpred2_Upper_99 = c(CI_99$R2_raw_EUR_RICE_vs_LDpred2_Upper,CI_99$R2_raw_SAS_RICE_vs_LDpred2_Upper,CI_99$R2_raw_AFR_RICE_vs_LDpred2_Upper,CI_99$R2_raw_AMR_RICE_vs_LDpred2_Upper),
                    R2_raw_RICE_vs_Lassosum2_Lower_99 = c(CI_99$R2_raw_EUR_RICE_vs_Lassosum2_Lower,CI_99$R2_raw_SAS_RICE_vs_Lassosum2_Lower,CI_99$R2_raw_AFR_RICE_vs_Lassosum2_Lower,CI_99$R2_raw_AMR_RICE_vs_Lassosum2_Lower),
                    R2_raw_RICE_vs_Lassosum2_Upper_99 = c(CI_99$R2_raw_EUR_RICE_vs_Lassosum2_Upper,CI_99$R2_raw_SAS_RICE_vs_Lassosum2_Upper,CI_99$R2_raw_AFR_RICE_vs_Lassosum2_Upper,CI_99$R2_raw_AMR_RICE_vs_Lassosum2_Upper),
                    R2_adjusted_RICE_vs_CT_Lower_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_CT_Lower,CI_99$R2_adjusted_SAS_RICE_vs_CT_Lower,CI_99$R2_adjusted_AFR_RICE_vs_CT_Lower,CI_99$R2_adjusted_AMR_RICE_vs_CT_Lower),
                    R2_adjusted_RICE_vs_CT_Upper_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_CT_Upper,CI_99$R2_adjusted_SAS_RICE_vs_CT_Upper,CI_99$R2_adjusted_AFR_RICE_vs_CT_Upper,CI_99$R2_adjusted_AMR_RICE_vs_CT_Upper),
                    R2_adjusted_RICE_vs_LDpred2_Lower_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_LDpred2_Lower,CI_99$R2_adjusted_SAS_RICE_vs_LDpred2_Lower,CI_99$R2_adjusted_AFR_RICE_vs_LDpred2_Lower,CI_99$R2_adjusted_AMR_RICE_vs_LDpred2_Lower),
                    R2_adjusted_RICE_vs_LDpred2_Upper_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_LDpred2_Upper,CI_99$R2_adjusted_SAS_RICE_vs_LDpred2_Upper,CI_99$R2_adjusted_AFR_RICE_vs_LDpred2_Upper,CI_99$R2_adjusted_AMR_RICE_vs_LDpred2_Upper),
                    R2_adjusted_RICE_vs_Lassosum2_Lower_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_Lassosum2_Lower,CI_99$R2_adjusted_SAS_RICE_vs_Lassosum2_Lower,CI_99$R2_adjusted_AFR_RICE_vs_Lassosum2_Lower,CI_99$R2_adjusted_AMR_RICE_vs_Lassosum2_Lower),
                    R2_adjusted_RICE_vs_Lassosum2_Upper_99 = c(CI_99$R2_adjusted_EUR_RICE_vs_Lassosum2_Upper,CI_99$R2_adjusted_SAS_RICE_vs_Lassosum2_Upper,CI_99$R2_adjusted_AFR_RICE_vs_Lassosum2_Upper,CI_99$R2_adjusted_AMR_RICE_vs_Lassosum2_Upper)) 
full_results <- left_join(full_results,CI_99)

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower_95 = full_results$beta_raw_Lower_95, upper_95 = full_results$beta_raw_Upper_95,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower_95 = full_results$beta_adjusted_Lower_95, upper_95 = full_results$beta_adjusted_Upper_95,method = full_results$Method,Standardization = "Using PCs 1-5"))

FigS14_Continuous <- full_results_stacked

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-CV"
full_results$group2 <- "RICE-CV"
full_results$p.signif_beta <- ""
full_results$p.signif_beta[full_results$Method == "RICE-CV"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV"] > 0,"**",""))
full_results$position <- NA
full_results$position[full_results$Method == "RICE-CV"] <- full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"] + 0.03
ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"],full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

FigS12_Continuous_Betas <- full_results

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV",]
full_results$Method[full_results$Method == "RICE-CV"] <- "RICE"
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

full_results$group1 <- "RICE"
full_results$group2 <- "RICE"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta2 <- ""

for(trait in c("BMI","Height","LDL","log(TG)","TC","HDL")){
  for(anc in c("AFR","EUR","SAS","AMR")){
    tmp <- full_results[full_results$ancestry == anc & full_results$trait == trait,]
    max_R2_notRICE <- max(tmp$R2_adjusted[tmp$Method != "RICE"])
    Best_Method <- tmp$Method[tmp$R2_adjusted == max_R2_notRICE]
    Improved_R2 <- round((tmp$R2_adjusted[tmp$Method == "RICE"]/max_R2_notRICE - 1)*100,digits = 2)
    
    if(Best_Method == "CT"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_99[tmp$Method == "RICE"] > 0,paste0("(***)"),ifelse(tmp$R2_adjusted_RICE_vs_CT_Lower_95[tmp$Method == "RICE"] > 0,paste0("(**)"),""))
    }else if(Best_Method == "LDpred2"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_LDpred2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else{
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_Lassosum2_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }
  }
}

full_results$position1 <- NA
full_results$position2 <- NA
full_results$position1[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.01
full_results$position2[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.07
ylim <- max(c(full_results$R2_adjusted)) + 0.08


full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

FigS12_Continuous_R2 <- full_results

#################################################################
### Supplementary Figure 13
#################################################################

FigS13_EUR <- NULL
FigS13_AMR <- NULL
FigS13_AFR <- NULL
FigS13_SAS <- NULL

for(trait in c("BMI","HDL","LDL","logTG","TC","Height")){
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
  pheno_validation$y_validation <- NA
  pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
  
  CV_RV_PRS_raw <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for(i in c("RV_PRS","CV_PRS")){
    tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      CV_RV_PRS_adjusted[,i] <- 0
    }else{
      CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
    }
  }
  
  
  CV_RV_PRS_adjusted$Common_Bin <- 9
  CV_RV_PRS_adjusted$Rare_Bin <- 9
  
  Common_quants <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(0,.1,.2,.3,.4,.6,.7,.8,.9))
  
  for(i in 1:8){
    CV_RV_PRS_adjusted$Common_Bin[CV_RV_PRS_adjusted$CV_PRS < unname(Common_quants[i + 1]) & CV_RV_PRS_adjusted$CV_PRS >= unname(Common_quants[i])] <- i
  }
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_EUR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_EUR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_EUR$y_validation <- scale(CV_RV_PRS_adjusted_EUR$y_validation)
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted_EUR[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin <- factor(CV_RV_PRS_adjusted_EUR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_AMR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AMR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_AMR$y_validation <- scale(CV_RV_PRS_adjusted_AMR$y_validation)
  
  CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted_AMR[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_AMR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AMR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_AFR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AFR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_AFR$y_validation <- scale(CV_RV_PRS_adjusted_AFR$y_validation)
  
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted_AFR[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_AFR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AFR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_SAS$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_SAS$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_SAS$y_validation <- scale(CV_RV_PRS_adjusted_SAS$y_validation)
  
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted_SAS[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_SAS$Rare_Bin <- factor(CV_RV_PRS_adjusted_SAS$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  CV_RV_PRS_adjusted_EUR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_EUR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_EUR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_EUR,mean)
  
  colnames(CV_RV_PRS_adjusted_EUR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_EUR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_EUR <- inner_join(CV_RV_PRS_adjusted_EUR,CV_RV_PRS_adjusted_EUR_se)
  
  CV_RV_PRS_adjusted_AMR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AMR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_AMR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AMR,mean)
  
  colnames(CV_RV_PRS_adjusted_AMR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_AMR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_AMR <- inner_join(CV_RV_PRS_adjusted_AMR,CV_RV_PRS_adjusted_AMR_se)
  
  CV_RV_PRS_adjusted_AFR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AFR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_AFR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AFR,mean)
  
  colnames(CV_RV_PRS_adjusted_AFR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_AFR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_AFR <- inner_join(CV_RV_PRS_adjusted_AFR,CV_RV_PRS_adjusted_AFR_se)
  
  CV_RV_PRS_adjusted_SAS_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_SAS,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_SAS <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_SAS,mean)
  
  colnames(CV_RV_PRS_adjusted_SAS) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_SAS_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_SAS <- inner_join(CV_RV_PRS_adjusted_SAS,CV_RV_PRS_adjusted_SAS_se)
  
  
  colnames(CV_RV_PRS_adjusted_EUR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_AMR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_AFR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_SAS) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  
  ymin_EUR <- round(min(c(CV_RV_PRS_adjusted_EUR$Mean - CV_RV_PRS_adjusted_EUR$SE)) - 0.05,2)
  ymax_EUR <- round(max(c(CV_RV_PRS_adjusted_EUR$Mean + CV_RV_PRS_adjusted_EUR$SE)) + 0.05,2)
  
  ymin_AMR <- round(min(c(CV_RV_PRS_adjusted_AMR$Mean - CV_RV_PRS_adjusted_AMR$SE)) - 0.05,2)
  ymax_AMR <- round(max(c(CV_RV_PRS_adjusted_AMR$Mean + CV_RV_PRS_adjusted_AMR$SE)) + 0.05,2)
  
  ymin_AFR <- round(min(c(CV_RV_PRS_adjusted_AFR$Mean - CV_RV_PRS_adjusted_AFR$SE)) - 0.05,2)
  ymax_AFR <- round(max(c(CV_RV_PRS_adjusted_AFR$Mean + CV_RV_PRS_adjusted_AFR$SE)) + 0.05,2)
  
  ymin_SAS <- round(min(c(CV_RV_PRS_adjusted_SAS$Mean - CV_RV_PRS_adjusted_SAS$SE)) - 0.05,2)
  ymax_SAS <- round(max(c(CV_RV_PRS_adjusted_SAS$Mean + CV_RV_PRS_adjusted_SAS$SE)) + 0.05,2)
  
  FigS13_EUR <- rbind(FigS13_EUR,cbind(CV_RV_PRS_adjusted_EUR,trait))
  FigS13_AMR <- rbind(FigS13_AMR,cbind(CV_RV_PRS_adjusted_AMR,trait))
  FigS13_AFR <- rbind(FigS13_AFR,cbind(CV_RV_PRS_adjusted_AFR,trait))
  FigS13_SAS <- rbind(FigS13_SAS,cbind(CV_RV_PRS_adjusted_SAS,trait))
}

#################################################################
### Figure 6 + Supplementary Figure 15
#################################################################


traits_cont <- c("BMI","TC","HDL","LDL","logTG","Height")
traits_bin  <- c("Asthma","Breast","CAD","Prostate","T2D")

## ------------------------------------------------------------
## Read Imputed + WES results (Continuous)
## ------------------------------------------------------------
Imputed_Results_Continuous <- NULL
for(trait in traits_cont){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  Imputed_Results_Continuous <- rbind(Imputed_Results_Continuous,
                                      rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
Imputed_Results_Continuous$Data_Type <- "Imputed + WES"

## ------------------------------------------------------------
## Read Imputed + WES results (Binary)
## ------------------------------------------------------------
Imputed_Results_Binary <- NULL
for(trait in traits_bin){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  Imputed_Results_Binary <- rbind(Imputed_Results_Binary,
                                  rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
Imputed_Results_Binary$Data_Type <- "Imputed + WES"

## ------------------------------------------------------------
## Read WGS results (Binary)
## ------------------------------------------------------------
WGS_Results_Binary <- NULL
for(trait in traits_bin){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  WGS_Results_Binary <- rbind(WGS_Results_Binary,
                              rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
WGS_Results_Binary$Data_Type <- "WGS"

## ------------------------------------------------------------
## Read WGS results (Continuous)
## ------------------------------------------------------------
WGS_Results_Continuous <- NULL
for(trait in traits_cont){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  WGS_Results_Continuous <- rbind(WGS_Results_Continuous,
                                  rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
WGS_Results_Continuous$Data_Type <- "WGS"

## ------------------------------------------------------------
## Combine + filter (DROP WES here)
## ------------------------------------------------------------
full_results_Continuous <- rbind(Imputed_Results_Continuous, WGS_Results_Continuous)
full_results_Continuous$Method[full_results_Continuous$Method == "CV"] <- "RICE-CV"
full_results_Continuous$Method[full_results_Continuous$Method == "RV"] <- "RICE-RV"
full_results_Continuous <- full_results_Continuous[full_results_Continuous$Method %in% c("RICE-CV","RICE-RV"),]
full_results_Continuous$Method_DataSource <- paste0(full_results_Continuous$Method," & ",full_results_Continuous$Data_Type)

full_results_Continuous <- full_results_Continuous[full_results_Continuous$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_Continuous$trait[full_results_Continuous$trait == "logTG"] <- "log(TG)"
full_results_Continuous$trait <- factor(full_results_Continuous$trait, levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))

full_results_Binary <- rbind(Imputed_Results_Binary, WGS_Results_Binary)
full_results_Binary$Method[full_results_Binary$Method == "CV"] <- "RICE-CV"
full_results_Binary$Method[full_results_Binary$Method == "RV"] <- "RICE-RV"
full_results_Binary <- full_results_Binary[full_results_Binary$Method %in% c("RICE-CV","RICE-RV"),]
full_results_Binary$Method_DataSource <- paste0(full_results_Binary$Method," & ",full_results_Binary$Data_Type)

full_results_Binary <- full_results_Binary[full_results_Binary$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_Binary$trait <- factor(full_results_Binary$trait, levels = c("Asthma","Breast","CAD","Prostate","T2D"))

## keep only these 4 legend entries
mds_levels_main <- c(
  "RICE-CV & Imputed + WES",
  "RICE-CV & WGS",
  "RICE-RV & Imputed + WES",
  "RICE-RV & WGS"
)
full_results_Continuous$Method_DataSource <- factor(full_results_Continuous$Method_DataSource, levels = mds_levels_main)
full_results_Binary$Method_DataSource <- factor(full_results_Binary$Method_DataSource, levels = mds_levels_main)

## clip negatives
full_results_Continuous$beta_adjusted[full_results_Continuous$beta_adjusted < 0] <- 0
full_results_Binary$beta_adjusted[full_results_Binary$beta_adjusted < 0] <- 0

## legend extraction ordering (alternate within dataset)
full_results_tmp <- full_results_Continuous
full_results_tmp$Method_DataSource <- factor(
  full_results_tmp$Method_DataSource,
  levels = c(
    "RICE-CV & Imputed + WES",
    "RICE-RV & Imputed + WES",
    "RICE-CV & WGS",
    "RICE-RV & WGS"
  )
)

ylim_continuous <- max(full_results_Continuous$beta_adjusted, na.rm = TRUE) + 0.03
ylim_binary <- max(full_results_Binary$beta_adjusted, na.rm = TRUE) + 0.03

Fig6_Continuous <- full_results_Continuous[full_results_Continuous$ancestry == "EUR",]
Fig6_Binary <- full_results_Binary[full_results_Binary$ancestry == "EUR",]

FigS15 <- full_results_Continuous[full_results_Continuous$ancestry != "EUR",]

#################################################################
### Supplementary Figure 16
#################################################################

full_results_binary <- NULL

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  Coding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/",trait,"_Coding_Best_Betas.csv"))
  Coding_Results$Method <- "Coding"
  Noncoding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/",trait,"_Noncoding_Best_Betas.csv"))
  Noncoding_Results$Method <- "Noncoding"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  full_results_binary <- rbind(full_results_binary,rbind(Coding_Results,Noncoding_Results,RICE_CV_Results,RICE_RV_Results))
}

full_results_binary <- full_results_binary[full_results_binary$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_binary <- full_results_binary[full_results_binary$Method %in% c("Coding","Noncoding","RICE-RV","RICE-CV"),]

full_results_binary$Method <- factor(full_results_binary$Method,levels = c("Coding","Noncoding","RICE-RV","RICE-CV"))
full_results_binary$trait <- factor(full_results_binary$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))
full_results_binary$ancestry <- factor(full_results_binary$ancestry,levels = c("AFR","AMR","EUR","SAS"))

full_results_binary$beta_adjusted[full_results_binary$beta_adjusted < 0] <- 0
full_results_binary$beta_raw[full_results_binary$beta_raw < 0] <- 0

ylim <- max(full_results_binary$beta_adjusted) + 0.05

FigS16_Binary <- full_results_binary

full_results_continuous <- NULL

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  Coding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"_Coding_Best_Betas.csv"))
  Coding_Results$Method <- "Coding"
  Noncoding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"_Noncoding_Best_Betas.csv"))
  Noncoding_Results$Method <- "Noncoding"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  full_results_continuous <- rbind(full_results_continuous,rbind(Coding_Results,Noncoding_Results,RICE_CV_Results,RICE_RV_Results))
}

full_results_continuous <- full_results_continuous[full_results_continuous$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_continuous <- full_results_continuous[full_results_continuous$Method %in% c("Coding","Noncoding","RICE-RV","RICE-CV"),]

full_results_continuous$Method <- factor(full_results_continuous$Method,levels = c("Coding","Noncoding","RICE-RV","RICE-CV"))
full_results_continuous$trait[full_results_continuous$trait == "logTG"] <- "log(TG)"
full_results_continuous$trait <- factor(full_results_continuous$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_continuous$ancestry <- factor(full_results_continuous$ancestry,levels = c("AFR","AMR","EUR","SAS"))


full_results_continuous$beta_adjusted[full_results_continuous$beta_adjusted < 0] <- 0
full_results_continuous$beta_raw[full_results_continuous$beta_raw < 0] <- 0

ylim <- max(full_results_continuous$beta_adjusted) + 0.05

FigS16_Continuous <- full_results_continuous

#################################################################
### Figure 7 + Supplementary Figure 19
#################################################################

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

Comparison_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Comparison_Boot.csv")

full_results <- rbind(CTSLEB_Results,PROSPER_Results,JointPRS_Results,RICE_CV_Results,RICE_RV_Results)
full_results_Boot <- rbind(CTSLEB_Boot_Results,PROSPER_Boot_Results,JointPRS_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results)

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
Comparison_Boot_Results$trait[Comparison_Boot_Results$trait == "logTG"] <- "log(TG)"
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

lower_95 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1)] <- paste0(colnames(lower_95)[-c(1)],"_Lower")
upper_95 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1)] <- paste0(colnames(upper_95)[-c(1)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_95)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CTSLEB_Lower_95 = c(CI_95$R2_raw_AFR_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_AMR_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_EAS_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_EUR_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_MID_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_SAS_RICE_vs_CTSLEB_Lower),
                    R2_raw_RICE_vs_CTSLEB_Upper_95 = c(CI_95$R2_raw_AFR_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_AMR_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_EAS_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_EUR_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_MID_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_SAS_RICE_vs_CTSLEB_Upper),
                    R2_raw_RICE_vs_PROSPER_Lower_95 = c(CI_95$R2_raw_AFR_RICE_vs_PROSPER_Lower,CI_95$R2_raw_AMR_RICE_vs_PROSPER_Lower,CI_95$R2_raw_EAS_RICE_vs_PROSPER_Lower,CI_95$R2_raw_EUR_RICE_vs_PROSPER_Lower,CI_95$R2_raw_MID_RICE_vs_PROSPER_Lower,CI_95$R2_raw_SAS_RICE_vs_PROSPER_Lower),
                    R2_raw_RICE_vs_PROSPER_Upper_95 = c(CI_95$R2_raw_AFR_RICE_vs_PROSPER_Upper,CI_95$R2_raw_AMR_RICE_vs_PROSPER_Upper,CI_95$R2_raw_EAS_RICE_vs_PROSPER_Upper,CI_95$R2_raw_EUR_RICE_vs_PROSPER_Upper,CI_95$R2_raw_MID_RICE_vs_PROSPER_Upper,CI_95$R2_raw_SAS_RICE_vs_PROSPER_Upper),
                    R2_raw_RICE_vs_JointPRS_Lower_95 = c(CI_95$R2_raw_AFR_RICE_vs_JointPRS_Lower,CI_95$R2_raw_AMR_RICE_vs_JointPRS_Lower,CI_95$R2_raw_EAS_RICE_vs_JointPRS_Lower,CI_95$R2_raw_EUR_RICE_vs_JointPRS_Lower,CI_95$R2_raw_MID_RICE_vs_JointPRS_Lower,CI_95$R2_raw_SAS_RICE_vs_JointPRS_Lower),
                    R2_raw_RICE_vs_JointPRS_Upper_95 = c(CI_95$R2_raw_AFR_RICE_vs_JointPRS_Upper,CI_95$R2_raw_AMR_RICE_vs_JointPRS_Upper,CI_95$R2_raw_EAS_RICE_vs_JointPRS_Upper,CI_95$R2_raw_EUR_RICE_vs_JointPRS_Upper,CI_95$R2_raw_MID_RICE_vs_JointPRS_Upper,CI_95$R2_raw_SAS_RICE_vs_JointPRS_Upper),
                    R2_adjusted_RICE_vs_CTSLEB_Lower_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_AMR_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_EAS_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_EUR_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_MID_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_SAS_RICE_vs_CTSLEB_Lower),
                    R2_adjusted_RICE_vs_CTSLEB_Upper_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_AMR_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_EAS_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_EUR_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_MID_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_SAS_RICE_vs_CTSLEB_Upper),
                    R2_adjusted_RICE_vs_PROSPER_Lower_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_AMR_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_EAS_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_EUR_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_MID_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_SAS_RICE_vs_PROSPER_Lower),
                    R2_adjusted_RICE_vs_PROSPER_Upper_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_AMR_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_EAS_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_EUR_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_MID_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_SAS_RICE_vs_PROSPER_Upper),
                    R2_adjusted_RICE_vs_JointPRS_Lower_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_AMR_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_EAS_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_EUR_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_MID_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_SAS_RICE_vs_JointPRS_Lower),
                    R2_adjusted_RICE_vs_JointPRS_Upper_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_AMR_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_EAS_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_EUR_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_MID_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_SAS_RICE_vs_JointPRS_Upper))
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_raw_Lower_99 = c(CI_99$beta_raw_AFR_boot_Lower,CI_99$beta_raw_AMR_boot_Lower,CI_99$beta_raw_EAS_boot_Lower,CI_99$beta_raw_EUR_boot_Lower,CI_99$beta_raw_MID_boot_Lower,CI_99$beta_raw_SAS_boot_Lower),
                    beta_raw_Upper_99 = c(CI_99$beta_raw_AFR_boot_Upper,CI_99$beta_raw_AMR_boot_Upper,CI_99$beta_raw_EAS_boot_Upper,CI_99$beta_raw_EUR_boot_Upper,CI_99$beta_raw_MID_boot_Upper,CI_99$beta_raw_SAS_boot_Upper),
                    R2_raw_Lower_99 = c(CI_99$R2_raw_AFR_boot_Lower,CI_99$R2_raw_AMR_boot_Lower,CI_99$R2_raw_EAS_boot_Lower,CI_99$R2_raw_EUR_boot_Lower,CI_99$R2_raw_MID_boot_Lower,CI_99$R2_raw_SAS_boot_Lower),
                    R2_raw_Upper_99 = c(CI_99$R2_raw_AFR_boot_Upper,CI_99$R2_raw_AMR_boot_Upper,CI_99$R2_raw_EAS_boot_Upper,CI_99$R2_raw_EUR_boot_Upper,CI_99$R2_raw_MID_boot_Upper,CI_99$R2_raw_SAS_boot_Upper),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower,CI_99$beta_adjusted_EAS_boot_Lower,CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_MID_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper,CI_99$beta_adjusted_EAS_boot_Upper,CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_MID_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper),
                    R2_adjusted_Lower_99 = c(CI_99$R2_adjusted_AFR_boot_Lower,CI_99$R2_adjusted_AMR_boot_Lower,CI_99$R2_adjusted_EAS_boot_Lower,CI_99$R2_adjusted_EUR_boot_Lower,CI_99$R2_adjusted_MID_boot_Lower,CI_99$R2_adjusted_SAS_boot_Lower),
                    R2_adjusted_Upper_99 = c(CI_99$R2_adjusted_AFR_boot_Upper,CI_99$R2_adjusted_AMR_boot_Upper,CI_99$R2_adjusted_EAS_boot_Upper,CI_99$R2_adjusted_EUR_boot_Upper,CI_99$R2_adjusted_MID_boot_Upper,CI_99$R2_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

lower_99 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.025)})
colnames(lower_99)[-c(1)] <- paste0(colnames(lower_99)[-c(1)],"_Lower")
upper_99 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.975)})
colnames(upper_99)[-c(1)] <- paste0(colnames(upper_99)[-c(1)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_99)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CTSLEB_Lower_99 = c(CI_99$R2_raw_AFR_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_AMR_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_EAS_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_EUR_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_MID_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_SAS_RICE_vs_CTSLEB_Lower),
                    R2_raw_RICE_vs_CTSLEB_Upper_99 = c(CI_99$R2_raw_AFR_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_AMR_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_EAS_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_EUR_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_MID_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_SAS_RICE_vs_CTSLEB_Upper),
                    R2_raw_RICE_vs_PROSPER_Lower_99 = c(CI_99$R2_raw_AFR_RICE_vs_PROSPER_Lower,CI_99$R2_raw_AMR_RICE_vs_PROSPER_Lower,CI_99$R2_raw_EAS_RICE_vs_PROSPER_Lower,CI_99$R2_raw_EUR_RICE_vs_PROSPER_Lower,CI_99$R2_raw_MID_RICE_vs_PROSPER_Lower,CI_99$R2_raw_SAS_RICE_vs_PROSPER_Lower),
                    R2_raw_RICE_vs_PROSPER_Upper_99 = c(CI_99$R2_raw_AFR_RICE_vs_PROSPER_Upper,CI_99$R2_raw_AMR_RICE_vs_PROSPER_Upper,CI_99$R2_raw_EAS_RICE_vs_PROSPER_Upper,CI_99$R2_raw_EUR_RICE_vs_PROSPER_Upper,CI_99$R2_raw_MID_RICE_vs_PROSPER_Upper,CI_99$R2_raw_SAS_RICE_vs_PROSPER_Upper),
                    R2_raw_RICE_vs_JointPRS_Lower_99 = c(CI_99$R2_raw_AFR_RICE_vs_JointPRS_Lower,CI_99$R2_raw_AMR_RICE_vs_JointPRS_Lower,CI_99$R2_raw_EAS_RICE_vs_JointPRS_Lower,CI_99$R2_raw_EUR_RICE_vs_JointPRS_Lower,CI_99$R2_raw_MID_RICE_vs_JointPRS_Lower,CI_99$R2_raw_SAS_RICE_vs_JointPRS_Lower),
                    R2_raw_RICE_vs_JointPRS_Upper_99 = c(CI_99$R2_raw_AFR_RICE_vs_JointPRS_Upper,CI_99$R2_raw_AMR_RICE_vs_JointPRS_Upper,CI_99$R2_raw_EAS_RICE_vs_JointPRS_Upper,CI_99$R2_raw_EUR_RICE_vs_JointPRS_Upper,CI_99$R2_raw_MID_RICE_vs_JointPRS_Upper,CI_99$R2_raw_SAS_RICE_vs_JointPRS_Upper),
                    R2_adjusted_RICE_vs_CTSLEB_Lower_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_AMR_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_EAS_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_EUR_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_MID_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_SAS_RICE_vs_CTSLEB_Lower),
                    R2_adjusted_RICE_vs_CTSLEB_Upper_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_AMR_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_EAS_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_EUR_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_MID_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_SAS_RICE_vs_CTSLEB_Upper),
                    R2_adjusted_RICE_vs_PROSPER_Lower_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_AMR_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_EAS_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_EUR_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_MID_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_SAS_RICE_vs_PROSPER_Lower),
                    R2_adjusted_RICE_vs_PROSPER_Upper_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_AMR_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_EAS_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_EUR_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_MID_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_SAS_RICE_vs_PROSPER_Upper),
                    R2_adjusted_RICE_vs_JointPRS_Lower_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_AMR_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_EAS_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_EUR_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_MID_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_SAS_RICE_vs_JointPRS_Lower),
                    R2_adjusted_RICE_vs_JointPRS_Upper_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_AMR_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_EAS_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_EUR_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_MID_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_SAS_RICE_vs_JointPRS_Upper))
full_results <- left_join(full_results,CI_99)

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower_95 = full_results$beta_raw_Lower_95, upper_95 = full_results$beta_raw_Upper_95,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower_95 = full_results$beta_adjusted_Lower_95, upper_95 = full_results$beta_adjusted_Upper_95,method = full_results$Method,Standardization = "Using PCs 1-5"))


full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-CV"
full_results$group2 <- "RICE-CV"
full_results$p.signif_beta <- ""
full_results$p.signif_beta[full_results$Method == "RICE-CV"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV"] > 0,"**",""))
full_results$position <- NA
full_results$position[full_results$Method == "RICE-CV"] <- full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"] + 0.03
ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"],full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-RV","RICE-CV"))

Fig7 <- full_results

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV",]
full_results$Method[full_results$Method == "RICE-CV"] <- "RICE"
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE"))

full_results$group1 <- "RICE"
full_results$group2 <- "RICE"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta2 <- ""

for(trait in c("BMI","Height","LDL","log(TG)","TC","HDL")){
  for(anc in c("AFR","EUR","SAS","AMR","MID","EAS")){
    tmp <- full_results[full_results$ancestry == anc & full_results$trait == trait,]
    max_R2_notRICE <- max(tmp$R2_adjusted[tmp$Method != "RICE"])
    Best_Method <- tmp$Method[tmp$R2_adjusted == max_R2_notRICE]
    Improved_R2 <- round((tmp$R2_adjusted[tmp$Method == "RICE"]/max_R2_notRICE - 1)*100,digits = 2)
    
    if(Best_Method == "CT-SLEB"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else if(Best_Method == "JointPRS"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else{
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }
  }
}

full_results$position1 <- NA
full_results$position2 <- NA
full_results$position1[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.01
full_results$position2[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.07
ylim <- max(c(full_results$R2_adjusted)) + 0.08

full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE"))

FigS19 <- full_results

#################################################################
### Figure 8 + Supplementary Figure 21
#################################################################

RICE_CV_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Results.csv")
RICE_CV_Results$Method <- "RICE-CV; Train (AoU) -> Validate (AoU)"
RICE_CV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Boot.csv")
RICE_CV_Boot_Results$Method <- "RICE-CV; Train (AoU) -> Validate (AoU)"
RICE_CV_Boot_Results <- RICE_CV_Boot_Results[,c("trait","beta_CV_adjusted_EUR_boot","beta_CV_adjusted_SAS_boot","beta_CV_adjusted_AMR_boot","beta_CV_adjusted_AFR_boot","Method")]
colnames(RICE_CV_Boot_Results) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")
RICE_RV_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Results.csv")
RICE_RV_Results$Method <- "RICE-RV; Train (AoU) -> Validate (AoU)"
RICE_RV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Boot.csv")
RICE_RV_Boot_Results$Method <- "RICE-RV; Train (AoU) -> Validate (AoU)"
RICE_RV_Boot_Results <- RICE_RV_Boot_Results[,c("trait","beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot","Method")]
colnames(RICE_RV_Boot_Results) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"Best_Betas_RICECV.csv"))
  tmp$Method <- "RICE-CV; Train (AoU) -> Validate (UKB)"
  RICE_CV_Results <- rbind(RICE_CV_Results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Bootstraps_RICECV.csv"))
  tmp$Method <- "RICE-CV; Train (AoU) -> Validate (UKB)"
  tmp <- tmp[,c("trait","beta_CV_adjusted_EUR_boot","beta_CV_adjusted_SAS_boot","beta_CV_adjusted_AMR_boot","beta_CV_adjusted_AFR_boot","Method")]
  colnames(tmp) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")
  RICE_CV_Boot_Results <- rbind(RICE_CV_Boot_Results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"Best_Betas_RICERV.csv"))
  tmp$Method <- "RICE-RV; Train (AoU) -> Validate (UKB)"
  RICE_RV_Results <- rbind(RICE_RV_Results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Bootstraps_RICERV.csv"))
  tmp$Method <- "RICE-RV; Train (AoU) -> Validate (UKB)"
  tmp <- tmp[,c("trait","beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot","Method")]
  colnames(tmp) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")
  RICE_RV_Boot_Results <- rbind(RICE_RV_Boot_Results,tmp)
}

full_results <- rbind(RICE_CV_Results,RICE_RV_Results)
full_results_Boot <- rbind(RICE_CV_Boot_Results,RICE_RV_Boot_Results)

full_results <- full_results[full_results$ancestry %in% c("AFR","AMR","EUR","SAS"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("RICE-RV; Train (AoU) -> Validate (UKB)","RICE-CV; Train (AoU) -> Validate (UKB)","RICE-RV; Train (AoU) -> Validate (AoU)","RICE-CV; Train (AoU) -> Validate (AoU)"))
full_results$Method1[full_results$Method1 == "RICE-RV; Train (AoU) -> Validate (UKB)"] <- "RICE-CV; Train (AoU) -> Validate (UKB)"
full_results$Method1[full_results$Method1 == "RICE-RV; Train (AoU) -> Validate (AoU)"] <- "RICE-CV; Train (AoU) -> Validate (AoU)"
full_results$Method1 <- factor(full_results$Method1,levels = c("RICE-CV; Train (AoU) -> Validate (UKB)","RICE-CV; Train (AoU) -> Validate (AoU)"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("AFR","AMR","EUR","SAS"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower,CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper,CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("AFR","AMR","EUR","SAS"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower,CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper,CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

full_results$Method_Dodge <- ifelse(str_detect(full_results$Method,"RICE-CV"),"RICE-CV","RICE-RV")

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-RV; Train (AoU) -> Validate (UKB)"
full_results$group2 <- "RICE-RV; Train (AoU) -> Validate (UKB)"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta1[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] > 0,"**",""))
full_results$position1 <- NA
full_results$position1[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] <- full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] + 0.03

full_results$group3 <- "RICE-RV; Train (AoU) -> Validate (AoU)"
full_results$group4 <- "RICE-RV; Train (AoU) -> Validate (AoU)"
full_results$p.signif_beta2 <- ""
full_results$p.signif_beta2[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] > 0,"**",""))
full_results$position2 <- NA
full_results$position2[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] <- full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] + 0.03


ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV; Train (AoU) -> Validate (UKB)"] + full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"],
              full_results$beta_adjusted[full_results$Method == "RICE-CV; Train (AoU) -> Validate (AoU)"] + full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"],
              full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("RICE-CV; Train (AoU) -> Validate (AoU)","RICE-CV; Train (AoU) -> Validate (UKB)","RICE-RV; Train (AoU) -> Validate (AoU)","RICE-RV; Train (AoU) -> Validate (UKB)"))

Fig8 <- full_results

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV; Train (AoU) -> Validate (AoU)",]
full_results <- full_results[full_results$Method != "RICE-RV; Train (AoU) -> Validate (UKB)",]
full_results$Method[full_results$Method == "RICE-CV; Train (AoU) -> Validate (AoU)"] <- "RICE; Train (AoU) -> Validate (AoU)"
full_results$Method[full_results$Method == "RICE-CV; Train (AoU) -> Validate (UKB)"] <- "RICE; Train (AoU) -> Validate (UKB)"
full_results$Method <- factor(full_results$Method,levels = c("RICE; Train (AoU) -> Validate (AoU)","RICE; Train (AoU) -> Validate (UKB)"))

ylim <- max(full_results$R2_adjusted) + 0.03

FigS21 <- full_results

#################################################################
### Figure R2
#################################################################

# for(trait in c("BMI","HDL","Height","LDL","logTG","TC")){
#   dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_",trait,".regenie"), sep="")
#   colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
#   dat$P <- 10^(-1*dat$LOG10P)
#   
#   dat1 <- dat[,c("CHROM","ID","REF","POS","ALT","LOG10P")]
#   colnames(dat1) <- c("CHR","SNP","REF","BP","ALT","P_UKB_PCs")
#   
#   dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_apc_",trait,".regenie"), sep="")
#   colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
#   dat$P <- 10^(-1*dat$LOG10P)
#   
#   dat2 <- dat[,c("CHROM","ID","REF","POS","ALT","LOG10P")]
#   colnames(dat2) <- c("CHR","SNP","REF","BP","ALT","P_1000G_PCs")
#   
#   dat <- inner_join(dat1,dat2)
#   
#   plot(dat1$P_UKB_PCs,dat2$P_1000G_PCs,xlab = "-log10(P) UKB PCs",ylab = "-log10(P) 1000G PCs")
#   abline(0, 1, col="red",lwd=1)
#   title(main = trait)
#   
#   print(trait)
#   print(summary(dat1$P_UKB_PCs - dat2$P_1000G_PCs))
#   
# }

#################################################################
### Figure R3
#################################################################

full_results <- NULL

for(trait in c("Height","BMI","TC","HDL","LDL","logTG")){
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_apc_Best_Betas.csv"))
  tmp$Method <- "CT_1000G_PCs"
  tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  full_results <- rbind(full_results,tmp) 
  
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  tmp$Method <- "CT_UKB_PCs"
  tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  full_results <- rbind(full_results,tmp) 
}

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0

FigR3 <- full_results 

#################################################################
### Figure R4
#################################################################

safe_ttest <- function(idx1, idx2, y) {
  # idx1, idx2 are index vectors; y is the numeric response vector
  if (length(idx1) < 2L || length(idx2) < 2L) {
    return(NA_real_)   # not enough data for a valid t-test
  }
  out <- tryCatch(
    t.test(y[idx1], y[idx2],
           alternative = "two.sided",
           var.equal   = FALSE)$p.value,
    error = function(e) NA_real_   # if t.test fails for any other reason
  )
  out
}

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
NRI_Data_Continuous <- NULL

FigR4_PlotB <- NULL
FigR4_PValues_PlotB <- NULL
FigR4_PlotA <- NULL

for (trait in continuous_traits) {
  
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  CV_PRS_Tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"))
  colnames(CV_PRS_Tune) <- c("IID","CV_PRS")
  pheno_tune <- inner_join(pheno_tune,CV_PRS_Tune)
  RV_PRS_Tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Tune.csv"))
  colnames(RV_PRS_Tune) <- c("IID","RV_PRS")
  pheno_tune <- inner_join(pheno_tune,RV_PRS_Tune)
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  ## ==== Residualize trait on covariates (robust to NA in covariates) ====
  model.null <- lm(
    as.formula(paste0(trait," ~ age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")),
    data = pheno_tune
  )
  pheno_tune$y_tune <- NA_real_
  # EDIT: align residuals with the exact rows used in the model
  rows_used_tune <- as.numeric(rownames(model.null$model))
  pheno_tune$y_tune[rows_used_tune] <- resid(model.null)
  
  model.null <- lm(
    as.formula(paste0(trait," ~ age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10")),
    data = pheno_validation
  )
  pheno_validation$y_validation <- NA_real_
  # EDIT: same robust residual alignment for validation
  rows_used_val <- as.numeric(rownames(model.null$model))
  pheno_validation$y_validation[rows_used_val] <- resid(model.null)
  
  ## ==== RICE model on tuning set ====
  RICE_Model <- lm(y_tune ~ CV_PRS + RV_PRS, data = pheno_tune)
  print(coef(RICE_Model))
  pheno_validation$PRS <- predict(RICE_Model, pheno_validation)
  
  ## Only keep non-missing phenotype for validation
  pheno_validation <- pheno_validation[!is.na(pheno_validation[, trait]), ]
  
  ## ==== PC-adjustment of PRS & CV_PRS ====
  CV_RV_PRS_raw      <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for (i in c("PRS","CV_PRS")) {
    tmp <- data.frame(
      y = CV_RV_PRS_adjusted[, i],
      CV_RV_PRS_adjusted[, c("pc1","pc2","pc3","pc4","pc5")]
    )
    mod <- lm(y ~ ., data = tmp)
    R   <- mod$residuals
    
    tmp2 <- data.frame(
      y = R^2,
      CV_RV_PRS_adjusted[, c("pc1","pc2","pc3","pc4","pc5")]
    )
    mod2  <- lm(y ~ ., data = tmp2)
    y_hat <- predict(mod2, tmp2)
    
    if (any(y_hat < 0, na.rm = TRUE)) {
      mod2  <- lm(y ~ 1, data = tmp2)
      y_hat <- predict(mod2, tmp2)
    }
    if (sum(sqrt(y_hat), na.rm = TRUE) == 0) {
      CV_RV_PRS_adjusted[, i] <- 0
    } else {
      CV_RV_PRS_adjusted[, i] <- R / sqrt(y_hat)
    }
  }
  
  # EDIT: actually use the PC-adjusted PRS values downstream
  pheno_validation$PRS    <- CV_RV_PRS_adjusted$PRS
  pheno_validation$CV_PRS <- CV_RV_PRS_adjusted$CV_PRS
  
  ## Clear per-trait container
  NRI_Data_Continuous <- NULL
  
  for (risk in c(0.05, 0.1)) {
    
    ## ==== Define "truth" and risk groups ====
    truth_HighRisk <- which(
      pheno_validation$y_validation > quantile(pheno_validation$y_validation, 0.9)
    )
    
    # EDIT: use >= for "high" cutpoints so everyone is in exactly one of low/high
    cv_cut  <- quantile(pheno_validation$CV_PRS,  0.9)
    prs_cut <- quantile(pheno_validation$PRS,  1 - risk)
    
    CV_HighRisk_RV_HighRisk <- which(
      pheno_validation$CV_PRS >= cv_cut &
        pheno_validation$PRS    >= prs_cut
    )
    CV_HighRisk_RV_NotHighRisk <- which(
      pheno_validation$CV_PRS >= cv_cut &
        pheno_validation$PRS    <  prs_cut
    )
    CV_NotHighRisk_RV_HighRisk <- which(
      pheno_validation$CV_PRS <  cv_cut &
        pheno_validation$PRS    >= prs_cut
    )
    CV_NotHighRisk_RV_NotHighRisk <- which(
      pheno_validation$CV_PRS <  cv_cut &
        pheno_validation$PRS    <  prs_cut
    )
    
    tmp <- data.frame(
      trait = trait,
      risk  = paste0(100 * risk, "%"),
      
      Percent_A = 100 * length(CV_NotHighRisk_RV_NotHighRisk) / nrow(pheno_validation),
      Mean_A    = mean(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk]),
      SE_A      = sd(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk]) /
        sqrt(length(CV_NotHighRisk_RV_NotHighRisk)),
      Total_Capture_A = sum(CV_NotHighRisk_RV_NotHighRisk %in% truth_HighRisk),
      
      Percent_B = 100 * length(CV_HighRisk_RV_NotHighRisk) / nrow(pheno_validation),
      Mean_B    = mean(pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk]),
      SE_B      = sd(pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk]) /
        sqrt(length(CV_HighRisk_RV_NotHighRisk)),
      Total_Capture_B = sum(CV_HighRisk_RV_NotHighRisk %in% truth_HighRisk),
      
      Percent_C = 100 * length(CV_NotHighRisk_RV_HighRisk) / nrow(pheno_validation),
      Mean_C    = mean(pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk]),
      SE_C      = sd(pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk]) /
        sqrt(length(CV_NotHighRisk_RV_HighRisk)),
      Total_Capture_C = sum(CV_NotHighRisk_RV_HighRisk %in% truth_HighRisk),
      
      Percent_D = 100 * length(CV_HighRisk_RV_HighRisk) / nrow(pheno_validation),
      Mean_D    = mean(pheno_validation$y_validation[CV_HighRisk_RV_HighRisk]),
      SE_D      = sd(pheno_validation$y_validation[CV_HighRisk_RV_HighRisk]) /
        sqrt(length(CV_HighRisk_RV_HighRisk)),
      Total_Capture_D = sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk),
      
      Mean_BD = mean(
        pheno_validation$y_validation[
          c(CV_HighRisk_RV_HighRisk, CV_HighRisk_RV_NotHighRisk)
        ]
      )
    )
    
    tmp$C_minus_A  <- (tmp$Mean_C  - tmp$Mean_A) / sd(pheno_validation$y_validation)
    tmp$B_minus_A  <- (tmp$Mean_B  - tmp$Mean_A) / sd(pheno_validation$y_validation)
    tmp$BD_minus_A <- (tmp$Mean_BD - tmp$Mean_A) / sd(pheno_validation$y_validation)
    
    tmp$A_vs_B <- safe_ttest(
      CV_NotHighRisk_RV_NotHighRisk,
      CV_HighRisk_RV_NotHighRisk,
      pheno_validation$y_validation
    )
    tmp$A_vs_C <- safe_ttest(
      CV_NotHighRisk_RV_NotHighRisk,
      CV_NotHighRisk_RV_HighRisk,
      pheno_validation$y_validation
    )
    tmp$A_vs_D <- safe_ttest(
      CV_NotHighRisk_RV_NotHighRisk,
      CV_HighRisk_RV_HighRisk,
      pheno_validation$y_validation
    )
    
    NRI_Data_Continuous <- rbind(NRI_Data_Continuous, tmp)
  } 
  
  ## ---------- 1. Proportion (donut) plot for risk = 10% ----------
  
  plot_data <- data.frame(
    Method = rep(
      c("Low CV PRS, Low PRS",
        "High CV PRS, Low PRS",
        "Low CV PRS, High PRS",
        "High CV PRS, High PRS"),
      each = nrow(NRI_Data_Continuous)
    ),
    trait = c(NRI_Data_Continuous$trait,
              NRI_Data_Continuous$trait,
              NRI_Data_Continuous$trait,
              NRI_Data_Continuous$trait),
    risk  = c(NRI_Data_Continuous$risk,
              NRI_Data_Continuous$risk,
              NRI_Data_Continuous$risk,
              NRI_Data_Continuous$risk),
    value = c(NRI_Data_Continuous$Total_Capture_A,
              NRI_Data_Continuous$Total_Capture_B,
              NRI_Data_Continuous$Total_Capture_C,
              NRI_Data_Continuous$Total_Capture_D)
  )
  
  plot_data$Method <- factor(
    plot_data$Method,
    levels = c("Low CV PRS, Low PRS",
               "High CV PRS, Low PRS",
               "Low CV PRS, High PRS",
               "High CV PRS, High PRS")
  )
  
  plot_data_sub <- plot_data[plot_data$trait == trait & plot_data$risk == "10%", ] %>%
    arrange(Method) %>%
    mutate(
      fraction = value / sum(value),
      ymax     = cumsum(fraction),
      ymin     = c(0, head(ymax, -1)),
      label_pos = (ymin + ymax) / 2,
      label     = paste0(round(100 * value / sum(value), 1), "%")
    )
  
  plot_data_sub$Method <- factor(
    plot_data_sub$Method,
    levels = c("Low CV PRS, Low PRS",
               "High CV PRS, Low PRS",
               "Low CV PRS, High PRS",
               "High CV PRS, High PRS")
  )
  
  FigR4_PlotA <- rbind(FigR4_PlotA,cbind(plot_data_sub,trait))
  
  ## ---------- 2. Mean standardized phenotype plot for risk = 10% ----------
  
  plot_data <- data.frame(
    Method = rep(
      c("Low CV PRS, Low PRS",
        "High CV PRS, Low PRS",
        "Low CV PRS, High PRS",
        "High CV PRS, High PRS"),
      each = nrow(NRI_Data_Continuous)
    ),
    trait = c(NRI_Data_Continuous$trait,
              NRI_Data_Continuous$trait,
              NRI_Data_Continuous$trait,
              NRI_Data_Continuous$trait),
    risk  = c(NRI_Data_Continuous$risk,
              NRI_Data_Continuous$risk,
              NRI_Data_Continuous$risk,
              NRI_Data_Continuous$risk),
    mean  = c(NRI_Data_Continuous$Mean_A,
              NRI_Data_Continuous$Mean_B,
              NRI_Data_Continuous$Mean_C,
              NRI_Data_Continuous$Mean_D),
    se    = c(NRI_Data_Continuous$SE_A,
              NRI_Data_Continuous$SE_B,
              NRI_Data_Continuous$SE_C,
              NRI_Data_Continuous$SE_D)
  )
  
  plot_data$Method <- factor(
    plot_data$Method,
    levels = c("Low CV PRS, Low PRS",
               "High CV PRS, Low PRS",
               "Low CV PRS, High PRS",
               "High CV PRS, High PRS")
  )
  
  ## p-values (A vs B, A vs C, A vs D), with BH adjustment
  p_raw <- c(
    NRI_Data_Continuous$A_vs_B,
    NRI_Data_Continuous$A_vs_C,
    NRI_Data_Continuous$A_vs_D
  )
  # EDIT: adjust for multiple testing (Benjamini–Hochberg)
  p_adj <- p.adjust(p_raw, method = "BH")
  
  if (trait %in% c("LDL","BMI","TC")) {
    stat.test <- data.frame(
      group1 = "Low CV PRS, Low PRS",
      group2 = rep(
        c("High CV PRS, Low PRS",
          "Low CV PRS, High PRS",
          "High CV PRS, High PRS"),
        each = nrow(NRI_Data_Continuous)
      ),
      trait = c(NRI_Data_Continuous$trait,
                NRI_Data_Continuous$trait,
                NRI_Data_Continuous$trait),
      risk  = c(NRI_Data_Continuous$risk,
                NRI_Data_Continuous$risk,
                NRI_Data_Continuous$risk),
      p.adj = signif(p_adj, 3),
      y.position = as.vector(
        outer(
          apply(
            cbind(NRI_Data_Continuous$Mean_B,
                  NRI_Data_Continuous$Mean_C,
                  NRI_Data_Continuous$Mean_D),
            1, max
          ),
          c(1.3, 1.4, 1.5),
          "*"
        )
      ),
      p.adj.signif = ifelse(
        p_adj > 0.05, "",
        ifelse(p_adj < 0.01, "**", "*")
      )
    )
  } else {
    stat.test <- data.frame(
      group1 = "Low CV PRS, Low PRS",
      group2 = rep(
        c("High CV PRS, Low PRS",
          "Low CV PRS, High PRS",
          "High CV PRS, High PRS"),
        each = nrow(NRI_Data_Continuous)
      ),
      trait = c(NRI_Data_Continuous$trait,
                NRI_Data_Continuous$trait,
                NRI_Data_Continuous$trait),
      risk  = c(NRI_Data_Continuous$risk,
                NRI_Data_Continuous$risk,
                NRI_Data_Continuous$risk),
      p.adj = signif(p_adj, 3),
      y.position = as.vector(
        outer(
          apply(
            cbind(NRI_Data_Continuous$Mean_B,
                  NRI_Data_Continuous$Mean_C,
                  NRI_Data_Continuous$Mean_D),
            1, max
          ),
          c(1.1, 1.2, 1.3),
          "*"
        )
      ),
      p.adj.signif = ifelse(
        p_adj > 0.05, "",
        ifelse(p_adj < 0.01, "**", "*")
      )
    )
  }
  
  # Keep only this trait + risk = 10% and non-NA p-values
  stat.test.sub <- stat.test[
    stat.test$trait == trait & stat.test$risk == "10%" & !is.na(stat.test$p.adj),
  ]
  
  FigR4_PlotB <- rbind(FigR4_PlotB,plot_data[plot_data$trait == trait & plot_data$risk == "10%", ])
  FigR4_PValues_PlotB <- rbind(FigR4_PValues_PlotB,stat.test.sub)
}




# ====================== Nature Communications Source Data export ======================
# - One workbook
# - One sheet per figure that EXISTS in the R environment
# - Expanded figure labels everywhere:
#     Fig3  -> "Figure 3"
#     FigS3 -> "Supplementary Figure 3"
#     FigR3 -> "Reviewer Figure 3"
# - "Contents" sheet columns:
#     Sheet name | Short description | Notes (missing figures)
# - Descriptions come from a MANUAL mapping in R (no .docx parsing)
# - Row names are NOT written

OUT_XLSX <- "NatureComm_SourceData.xlsx"   # <- change output filename/path if desired

# ---- MANUAL mapping: edit these descriptions to match your manuscript/supp docs ----
# Keys MUST be base IDs only: "Fig1", "Fig2", "FigS1", "FigR1", etc.
FIG_DESCRIPTIONS <- c(
  # Main figures
  "Fig1" = "Overview of the RICE framework for polygenic risk prediction.",
  "Fig2" = "Comparison of standardized PRSs from RICE-CV and RICE-RV for high-density lipoproteins cholesterol (HDL).",
  "Fig3" = "Simulation results comparing the predictive performance of PRSs for four ancestral groups from the UK Biobank (UKB).",
  "Fig4" = "Predictive performance of ancestry-adjusted PRSs for continuous traits across four ancestral groups from UK Biobank (UKB) imputed + whole-exome sequencing (WES) data.",
  "Fig5" = "Joint stratification of common and rare variant PRS identifies discrepant high-risk individuals.",
  "Fig6" = "Comparison of ancestry-adjusted PRSs from RICE-CV and RICE-RV using UK Biobank (UKB) in two configurations: Imputed + WES (imputed genotypes for common variants and whole exome sequencing data for rare variants) and WGS (whole-genome sequencing data for both common and rare variants).",
  "Fig7" = "Predictive performance of ancestry-adjusted PRSs for continuous traits across six ancestral groups from the All of Us (AoU) whole-exome sequencing (WES) data.",
  "Fig8" = "Predictive performance of RICE trained on All of Us (AoU) data and evaluated on both AoU and UK Biobank (UKB) validation datasets.",
  
  # Supplementary figures
  "FigS1"  = "Simulation results comparing the predictive performance of PRSs for four ancestral groups from the UK Biobank (UKB).",
  "FigS2"  = "Average heritability of rare variant burden scores for each ancestry by simulation design.",
  "FigS3"  = "Manhattan and QQ plots for UKB Imputed + WES.",
  "FigS4"  = "Rare variant QQ plots for UKB Imputed + WES.",
  "FigS5"  = "Manhattan and QQ plots for UKB WGS.",
  "FigS6"  = "Rare variant QQ plots for UKB WGS.",
  "FigS7"  = "Predictive performance of ancestry-adjusted PRSs for 11 traits across four ancestral groups from UK Biobank (UKB) imputed + whole exome sequencing data (WES) data.",
  "FigS8"  = "Relationship between common and rare variant PRSs and standardized traits across four ancestral groups from UK Biobank (UKB) imputed + whole-exome sequencing (WES) data.",
  "FigS9"  = "Joint stratification of common and rare variant PRS for standardized traits across four ancestral groups from UK Biobank (UKB) imputed + whole-exome sequencing (WES) data.",
  "FigS10" = "Predictive performance of RICE-RV constructed using high-penetrance genes.",
  "FigS11" = "Predictive performance of RICE-RV constructed using different p-value thresholds.",
  "FigS12" = "Predictive performance of ancestry-adjusted PRSs for 11 traits across four ancestral groups from UK Biobank (UKB) whole genome sequencing data (WGS) data.",
  "FigS13" = "Relationship between common and rare variant PRSs and standardized traits across four ancestral groups from UK Biobank (UKB) whole-genome sequencing (WGS) data.",
  "FigS14" = "Predictive performance of PRSs standardized within genetically-inferred ancestries or using the first five principal components (Methods) for six continuous traits and five binary traits across four ancestral groups from UK Biobank (UKB) whole-genome sequencing (WGS) data.",
  "FigS15" = "Comparison of ancestry-adjusted PRSs from RICE-CV and RICE-RV using UK Biobank (UKB) whole-exome sequencing (WES), imputed + WES, and whole-genome sequencing (WGS) data for individuals of either African, Admixed American, or South Asian ancestry.",
  "FigS16" = "Comparison of ancestry-adjusted PRSs from RICE-CV, RICE-RV, RICE-RV constructed with only coding genes (Coding), and RICE-RV constructed using only noncoding genes (Noncoding) using UK Biobank (UKB) whole-genome sequencing (WGS) data for individuals of either African, Admixed American, or South Asian ancestry.",
  "FigS17" = "Manhattan and QQ plots for All of Us.",
  "FigS18" = "Rare variant QQ plots for All of Us.",
  "FigS19" = "Predictive performance of ancestry-adjusted PRSs for six continuous traits across six ancestral groups from All of Us (AoU) data.",
  "FigS20" = "Relationship between common and rare variant PRSs and standardized traits across six ancestral groups from All of Us (AoU) data.",
  "FigS21" = "Assessing the prediction performance of RICE trained on All of Us (AoU) and validated on UK Biobank Imputed + whole-exome sequencing (WES) data."
)

# Optional overrides to force the TOC to include up to these numbers even if absent
MAX_MAIN_OVERRIDE <- NA_integer_  # e.g., 8L  for Figure 1..8
MAX_SUPP_OVERRIDE <- NA_integer_  # e.g., 25L for Supplementary Figure 1..25
MAX_REV_OVERRIDE  <- NA_integer_  # e.g., 6L  for Reviewer Figure 1..6

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required. Install it with: install.packages('openxlsx')")
}
library(openxlsx)

sd_env <- .GlobalEnv

# --- helpers ------------------------------------------------------------------------
sd_extract_num <- function(fig_id) as.integer(sub("^Fig(S|R)?", "", fig_id))

sd_pretty_fig <- function(fig_id) {
  n <- sd_extract_num(fig_id)
  if (grepl("^FigS", fig_id)) return(paste0("Supplementary Figure ", n))
  if (grepl("^FigR", fig_id)) return(paste0("Reviewer Figure ", n))
  paste0("Figure ", n)
}

sd_desc_for <- function(fig_id) {
  if (!is.null(names(FIG_DESCRIPTIONS)) && fig_id %in% names(FIG_DESCRIPTIONS)) {
    return(unname(FIG_DESCRIPTIONS[[fig_id]]))
  }
  ""
}

sd_as_tables <- function(x) {
  if (is.data.frame(x)) return(list(data = x))
  if (is.matrix(x))     return(list(data = as.data.frame(x)))
  if (is.atomic(x) && is.null(dim(x))) {
    return(list(data = data.frame(value = x, stringsAsFactors = FALSE)))
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) > 0 && all(vapply(x, function(e) is.data.frame(e) || is.matrix(e), logical(1)))) {
      out <- lapply(x, function(e) if (is.matrix(e)) as.data.frame(e) else e)
      if (is.null(names(out)) || any(names(out) == "")) names(out) <- paste0("element_", seq_along(out))
      return(out)
    }
    df <- tryCatch(as.data.frame(x), error = function(e) NULL)
    if (!is.null(df)) return(list(data = df))
  }
  list(structure = data.frame(value = capture.output(str(x)), stringsAsFactors = FALSE))
}

# Warn if mapping keys look wrong
bad_keys <- names(FIG_DESCRIPTIONS)[!grepl("^Fig(S|R)?[0-9]+$", names(FIG_DESCRIPTIONS))]
if (length(bad_keys)) {
  warning("FIG_DESCRIPTIONS has invalid keys (must look like Fig3, FigS12, FigR4): ",
          paste(bad_keys, collapse = ", "))
}

# --- discover Fig* objects ----------------------------------------------------------
sd_obj_names <- ls(envir = sd_env, pattern = "^Fig", all.names = TRUE)
sd_obj_names <- sd_obj_names[grepl("^Fig(S|R)?[0-9]+", sd_obj_names)]

# drop functions/environments
sd_obj_names <- sd_obj_names[!vapply(sd_obj_names, function(nm) {
  x <- get(nm, envir = sd_env)
  is.function(x) || is.environment(x)
}, logical(1))]

# Base figure like Fig3 / FigS12 / FigR4
sd_base_fig <- sub("^(Fig(?:S|R)?[0-9]+).*", "\\1", sd_obj_names)
sd_by_fig <- split(sd_obj_names, sd_base_fig)

# --- decide which figures to list in TOC (align order: main, supp, reviewer) --------
sd_all_ids_from_obj <- names(sd_by_fig)
sd_all_ids_from_map <- names(FIG_DESCRIPTIONS)

sd_main_ids <- unique(c(sd_all_ids_from_obj, sd_all_ids_from_map))[grepl("^Fig[0-9]+$", unique(c(sd_all_ids_from_obj, sd_all_ids_from_map)))]
sd_supp_ids <- unique(c(sd_all_ids_from_obj, sd_all_ids_from_map))[grepl("^FigS[0-9]+$", unique(c(sd_all_ids_from_obj, sd_all_ids_from_map)))]
sd_rev_ids  <- unique(c(sd_all_ids_from_obj, sd_all_ids_from_map))[grepl("^FigR[0-9]+$", unique(c(sd_all_ids_from_obj, sd_all_ids_from_map)))]

sd_max_main <- if (length(sd_main_ids)) max(sd_extract_num(sd_main_ids), na.rm = TRUE) else 0L
sd_max_supp <- if (length(sd_supp_ids)) max(sd_extract_num(sd_supp_ids), na.rm = TRUE) else 0L
sd_max_rev  <- if (length(sd_rev_ids))  max(sd_extract_num(sd_rev_ids),  na.rm = TRUE) else 0L

if (!is.na(MAX_MAIN_OVERRIDE)) sd_max_main <- MAX_MAIN_OVERRIDE
if (!is.na(MAX_SUPP_OVERRIDE)) sd_max_supp <- MAX_SUPP_OVERRIDE
if (!is.na(MAX_REV_OVERRIDE))  sd_max_rev  <- MAX_REV_OVERRIDE

sd_all_fig_ids <- c(
  if (sd_max_main > 0) paste0("Fig",  seq_len(sd_max_main)) else character(0),
  if (sd_max_supp > 0) paste0("FigS", seq_len(sd_max_supp)) else character(0),
  if (sd_max_rev  > 0) paste0("FigR", seq_len(sd_max_rev))  else character(0)
)

sd_present <- sd_all_fig_ids %in% names(sd_by_fig)

# --- Contents (Sheet name | Short description | Notes) -------------------------------
sd_toc <- data.frame(
  `Sheet name`        = vapply(sd_all_fig_ids, sd_pretty_fig, character(1)),
  `Short description` = vapply(sd_all_fig_ids, sd_desc_for, character(1)),
  Notes               = ifelse(sd_present, "", "Missing (no Fig* source-data object found)"),
  stringsAsFactors = FALSE
)

# --- Build workbook -----------------------------------------------------------------
sd_wb <- createWorkbook()

addWorksheet(sd_wb, "Contents")
writeData(sd_wb, "Contents", sd_toc)
freezePane(sd_wb, "Contents", firstRow = TRUE)
setColWidths(sd_wb, "Contents", cols = 1:ncol(sd_toc), widths = "auto")

hdr_style   <- createStyle(textDecoration = "bold")
title_style <- createStyle(textDecoration = "bold", fontSize = 14)
subtle_style <- createStyle(fontSize = 10)

addStyle(sd_wb, "Contents", hdr_style, rows = 1, cols = 1:ncol(sd_toc), gridExpand = TRUE)

# Create sheets ONLY for figures that exist
for (sd_fig_id in sd_all_fig_ids[sd_present]) {
  sd_sheet <- sd_pretty_fig(sd_fig_id)
  addWorksheet(sd_wb, sd_sheet)
  
  # Top-of-sheet header
  writeData(sd_wb, sd_sheet, x = sd_sheet, startRow = 1, startCol = 1, colNames = FALSE)
  addStyle(sd_wb, sd_sheet, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  
  sd_desc <- sd_desc_for(sd_fig_id)
  if (nzchar(sd_desc)) {
    writeData(sd_wb, sd_sheet, x = sd_desc, startRow = 2, startCol = 1, colNames = FALSE)
    addStyle(sd_wb, sd_sheet, subtle_style, rows = 2, cols = 1, gridExpand = TRUE)
  }
  
  sd_r <- 4
  
  for (sd_obj_name in sd_by_fig[[sd_fig_id]]) {
    sd_obj <- get(sd_obj_name, envir = sd_env)
    sd_tables <- sd_as_tables(sd_obj)
    
    for (sd_tbl_name in names(sd_tables)) {
      sd_label <- if (sd_tbl_name == "data") {
        paste0(sd_sheet, " — ", sd_obj_name)
      } else {
        paste0(sd_sheet, " — ", sd_obj_name, " / ", sd_tbl_name)
      }
      
      writeData(sd_wb, sd_sheet, x = sd_label, startRow = sd_r, startCol = 1, colNames = FALSE)
      addStyle(sd_wb, sd_sheet, hdr_style, rows = sd_r, cols = 1, gridExpand = TRUE)
      sd_r <- sd_r + 1
      
      sd_df <- sd_tables[[sd_tbl_name]]
      writeData(sd_wb, sd_sheet, x = sd_df, startRow = sd_r, startCol = 1, rowNames = FALSE)
      sd_r <- sd_r + max(nrow(sd_df), 1) + 3
    }
  }
  
  setColWidths(sd_wb, sd_sheet, cols = 1:50, widths = "auto")
}

saveWorkbook(sd_wb, OUT_XLSX, overwrite = TRUE)
message("Wrote source data workbook: ", normalizePath(OUT_XLSX, winslash = "/", mustWork = FALSE))
# ================================================================================ 