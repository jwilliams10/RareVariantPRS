rm(list = ls())

library(ggplot2)
library(ggpubr)
library(dplyr)

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

results_70 <- NULL
results_70_CIs <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv"))
  
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
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_70 <- rbind(results_70,betas_tmp)
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_70_CIs <- rbind(results_70_CIs,CIs_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","i","Y_train","index_mat")))
}

results_70 <- inner_join(results_70,index_mat)
results_70$Causal_Prop <- as.character(results_70$Causal_Prop)
results_70$Causal_Prop[results_70$Causal_Prop == "5e-04"] <- "0.0005"
results_70$Causal_Prop <- paste0("Causal Prop. ",results_70$Causal_Prop)

results_70$Scale <- as.character(results_70$Scale)
results_70$Scale[results_70$Scale == "0"] <- "Unscaled"
results_70$Scale[results_70$Scale == "1"] <- "Scaled"

results_70 <- data.frame(Scale = results_70$Scale, Causal_Prop = results_70$Causal_Prop, Method = results_70$Method,Ancestry = results_70$ancestry,
                         Beta = results_70$beta_adjusted,SE_Beta = results_70$beta_se_adjusted,R2 = results_70$R2_adjusted,SE_R2 = results_70$R2_se_adjusted)
results_70$Train_Size <- nrow(Y_train[[1]])

results_70_CIs <- inner_join(results_70_CIs,index_mat)
results_70_CIs$Causal_Prop <- as.character(results_70_CIs$Causal_Prop)
results_70_CIs$Causal_Prop[results_70_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_70_CIs$Causal_Prop <- paste0("Causal Prop. ",results_70_CIs$Causal_Prop)

results_70_CIs$Scale <- as.character(results_70_CIs$Scale)
results_70_CIs$Scale[results_70_CIs$Scale == "0"] <- "Unscaled"
results_70_CIs$Scale[results_70_CIs$Scale == "1"] <- "Scaled"

results_70_CIs$Train_Size <- nrow(Y_train[[1]])







load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_35 <- NULL
results_35_CIs <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv"))

  
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
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_35 <- rbind(results_35,betas_tmp)
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_35_CIs <- rbind(results_35_CIs,CIs_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_35","results_35_CIs","i","Y_train","index_mat")))
}

results_35 <- inner_join(results_35,index_mat)
results_35$Causal_Prop <- as.character(results_35$Causal_Prop)
results_35$Causal_Prop[results_35$Causal_Prop == "5e-04"] <- "0.0005"
results_35$Causal_Prop <- paste0("Causal Prop. ",results_35$Causal_Prop)

results_35$Scale <- as.character(results_35$Scale)
results_35$Scale[results_35$Scale == "0"] <- "Unscaled"
results_35$Scale[results_35$Scale == "1"] <- "Scaled"

results_35 <- data.frame(Scale = results_35$Scale, Causal_Prop = results_35$Causal_Prop, Method = results_35$Method,Ancestry = results_35$ancestry,
                         Beta = results_35$beta_adjusted,SE_Beta = results_35$beta_se_adjusted,R2 = results_35$R2_adjusted,SE_R2 = results_35$R2_se_adjusted)
results_35$Train_Size <- nrow(Y_train[[1]])

results_35_CIs <- inner_join(results_35_CIs,index_mat)
results_35_CIs$Causal_Prop <- as.character(results_35_CIs$Causal_Prop)
results_35_CIs$Causal_Prop[results_35_CIs$Causal_Prop == "5e-04"] <- "0.0005"
results_35_CIs$Causal_Prop <- paste0("Causal Prop. ",results_35_CIs$Causal_Prop)

results_35_CIs$Scale <- as.character(results_35_CIs$Scale)
results_35_CIs$Scale[results_35_CIs$Scale == "0"] <- "Unscaled"
results_35_CIs$Scale[results_35_CIs$Scale == "1"] <- "Scaled"

results_35_CIs$Train_Size <- nrow(Y_train[[1]])





load("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_rareprop_70 <- NULL
results_rareprop_70_CI <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv"))

  
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
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_rareprop_70 <- rbind(results_rareprop_70,betas_tmp)
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_rareprop_70_CI <- rbind(results_rareprop_70_CI,CIs_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_35","results_35_CIs","results_rareprop_70","results_rareprop_70_CI","i","Y_train","index_mat")))
}

results_rareprop_70 <- inner_join(results_rareprop_70,index_mat)
results_rareprop_70$Causal_Prop <- as.character(results_rareprop_70$Causal_Prop)
results_rareprop_70$Causal_Prop[results_rareprop_70$Causal_Prop == "5e-04"] <- "0.0005"
results_rareprop_70$Causal_Prop <- paste0("Causal Prop. ",results_rareprop_70$Causal_Prop)

results_rareprop_70$Scale <- as.character(results_rareprop_70$Scale)
results_rareprop_70$Scale[results_rareprop_70$Scale == "0"] <- "Unscaled"
results_rareprop_70$Scale[results_rareprop_70$Scale == "1"] <- "Scaled"

results_rareprop_70 <- data.frame(Scale = results_rareprop_70$Scale, Causal_Prop = results_rareprop_70$Causal_Prop, Method = results_rareprop_70$Method,Ancestry = results_rareprop_70$ancestry,
                         Beta = results_rareprop_70$beta_adjusted,SE_Beta = results_rareprop_70$beta_se_adjusted,R2 = results_rareprop_70$R2_adjusted,SE_R2 = results_rareprop_70$R2_se_adjusted)
results_rareprop_70$Train_Size <- nrow(Y_train[[1]])

results_rareprop_70_CI <- inner_join(results_rareprop_70_CI,index_mat)
results_rareprop_70_CI$Causal_Prop <- as.character(results_rareprop_70_CI$Causal_Prop)
results_rareprop_70_CI$Causal_Prop[results_rareprop_70_CI$Causal_Prop == "5e-04"] <- "0.0005"
results_rareprop_70_CI$Causal_Prop <- paste0("Causal Prop. ",results_rareprop_70_CI$Causal_Prop)

results_rareprop_70_CI$Scale <- as.character(results_rareprop_70_CI$Scale)
results_rareprop_70_CI$Scale[results_rareprop_70_CI$Scale == "0"] <- "Unscaled"
results_rareprop_70_CI$Scale[results_rareprop_70_CI$Scale == "1"] <- "Scaled"

results_rareprop_70_CI$Train_Size <- nrow(Y_train[[1]])







load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_rareprop_35 <- NULL
results_rareprop_35_CI <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Bootstraps_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/RV_",i,"_Bootstraps.csv"))

  
  lower_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.025)})
  lower_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_95 = lower_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_95 = lower_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_95 = lower_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_95 = lower_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_95 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.975)})
  upper_95 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_95 = upper_95[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_95 = upper_95[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_95 = upper_95[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_95 = upper_95[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  lower_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.005)})
  lower_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),Beta_Raw_Lower_99 = lower_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Lower_99 = lower_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Lower_99 = lower_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Lower_99 = lower_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  upper_99 <- apply(Bootstraps_RICERV,2,function(x){quantile(x[!is.na(x) & x <= 1 & x >= -1],0.995)})
  upper_99 <- data.frame(i = Bootstraps_RICERV$i[1],Ancestry = c("EUR","SAS","AMR","AFR"),Beta_Raw_Upper_99 = upper_99[c("beta_RV_raw_EUR_boot","beta_RV_raw_SAS_boot","beta_RV_raw_AMR_boot","beta_RV_raw_AFR_boot")],
                         Beta_Adjusted_Upper_99 = upper_99[c("beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot")],
                         R2_Raw_Upper_99 = upper_99[c("R2_raw_EUR_boot","R2_raw_SAS_boot","R2_raw_AMR_boot","R2_raw_AFR_boot")],
                         R2_Adjusted_Upper_99 = upper_99[c("R2_adjusted_EUR_boot","R2_adjusted_SAS_boot","R2_adjusted_AMR_boot","R2_adjusted_AFR_boot")])
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_rareprop_35 <- rbind(results_rareprop_35,betas_tmp)
  
  CIs_tmp <- inner_join(lower_95,upper_95)
  CIs_tmp <- inner_join(CIs_tmp,lower_99)
  CIs_tmp <- inner_join(CIs_tmp,upper_99)
  results_rareprop_35_CI <- rbind(results_rareprop_35_CI,CIs_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_70_CIs","results_35","results_35_CIs","results_rareprop_70","results_rareprop_70_CI","results_rareprop_35","results_rareprop_35_CI","i","Y_train","index_mat")))
}

results_rareprop_35 <- inner_join(results_rareprop_35,index_mat)
results_rareprop_35$Causal_Prop <- as.character(results_rareprop_35$Causal_Prop)
results_rareprop_35$Causal_Prop[results_rareprop_35$Causal_Prop == "5e-04"] <- "0.0005"
results_rareprop_35$Causal_Prop <- paste0("Causal Prop. ",results_rareprop_35$Causal_Prop)

results_rareprop_35$Scale <- as.character(results_rareprop_35$Scale)
results_rareprop_35$Scale[results_rareprop_35$Scale == "0"] <- "Unscaled"
results_rareprop_35$Scale[results_rareprop_35$Scale == "1"] <- "Scaled"

results_rareprop_35 <- data.frame(Scale = results_rareprop_35$Scale, Causal_Prop = results_rareprop_35$Causal_Prop, Method = results_rareprop_35$Method,Ancestry = results_rareprop_35$ancestry,
                                  Beta = results_rareprop_35$beta_adjusted,SE_Beta = results_rareprop_35$beta_se_adjusted,R2 = results_rareprop_35$R2_adjusted,SE_R2 = results_rareprop_35$R2_se_adjusted)
results_rareprop_35$Train_Size <- nrow(Y_train[[1]])

results_rareprop_35_CI <- inner_join(results_rareprop_35_CI,index_mat)
results_rareprop_35_CI$Causal_Prop <- as.character(results_rareprop_35_CI$Causal_Prop)
results_rareprop_35_CI$Causal_Prop[results_rareprop_35_CI$Causal_Prop == "5e-04"] <- "0.0005"
results_rareprop_35_CI$Causal_Prop <- paste0("Causal Prop. ",results_rareprop_35_CI$Causal_Prop)

results_rareprop_35_CI$Scale <- as.character(results_rareprop_35_CI$Scale)
results_rareprop_35_CI$Scale[results_rareprop_35_CI$Scale == "0"] <- "Unscaled"
results_rareprop_35_CI$Scale[results_rareprop_35_CI$Scale == "1"] <- "Scaled"

results_rareprop_35_CI$Train_Size <- nrow(Y_train[[1]])







results <- rbind(results_35,results_70)
results_CI <- rbind(results_35_CIs,results_70_CIs)
results_rareprop <- rbind(results_rareprop_35,results_rareprop_70)
results_rareprop_CI <- rbind(results_rareprop_35_CI,results_rareprop_70_CI)

results$Train_Size <- format(results$Train_Size,big.mark=",", trim=TRUE)
results_rareprop$Train_Size <- format(results_rareprop$Train_Size,big.mark=",", trim=TRUE)

results_CI$Train_Size <- format(results_CI$Train_Size,big.mark=",", trim=TRUE)
results_rareprop_CI$Train_Size <- format(results_rareprop_CI$Train_Size,big.mark=",", trim=TRUE)

results$Train_Size <- paste0("n = ",results$Train_Size)
results_rareprop$Train_Size <- paste0("n = ",results_rareprop$Train_Size)

results_CI$Train_Size <- paste0("n = ",results_CI$Train_Size)
results_rareprop_CI$Train_Size <- paste0("n = ",results_rareprop_CI$Train_Size)

rm(list=setdiff(ls(), c("results","results_CI","results_rareprop","results_rareprop_CI")))

results$Beta[results$Method %in% c("LDPred","LASSOSum")] <- -1*results$Beta[results$Method %in% c("LDPred","LASSOSum")]
results$Beta[results$Beta < 0] <- 0
results$Beta[results$Beta > 1] <- 0

results_rareprop$Beta[results_rareprop$Method %in% c("LDPred","LASSOSum")] <- -1*results_rareprop$Beta[results_rareprop$Method %in% c("LDPred","LASSOSum")]
results_rareprop$Beta[results_rareprop$Beta < 0] <- 0
results_rareprop$Beta[results_rareprop$Beta > 1] <- 0

results <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results,mean)
results_rareprop <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results_rareprop,mean)

results_CI <- aggregate(.~Scale + Causal_Prop + Train_Size + Ancestry,data = subset(results_CI,select = -c(i)),mean)
results_rareprop_CI <- aggregate(.~Scale + Causal_Prop + Train_Size + Ancestry,data = subset(results_rareprop_CI,select = -c(i)),mean)

overall_results <- results[results$Method %in% c("CT","LASSOSum","LDPred","RICE-CV","RICE-RV"),]
overall_results_rareprop <- results_rareprop[results_rareprop$Method %in% c("CT","LASSOSum","LDPred","RICE-CV","RICE-RV"),]

overall_results <- overall_results[overall_results$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results$Method[overall_results$Method == "RICE-CV"] <- "RICE-CV" 
overall_results$Method[overall_results$Method == "RICE-RV"] <- "RICE-RV" 
overall_results$Method[overall_results$Method == "LDPred"] <- "LDpred2"
overall_results$Method[overall_results$Method == "LASSOSum"] <- "Lassosum2"

overall_results$Method1 <- overall_results$Method
overall_results$Method <- factor(overall_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
overall_results$Method1[overall_results$Method1 == "RICE-RV"] <- "RICE-CV"
overall_results$Method1 <- factor(overall_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))
overall_results$Ancestry <- factor(overall_results$Ancestry,levels = c("AFR","AMR","EUR","SAS"))

overall_results_rareprop <- overall_results_rareprop[overall_results_rareprop$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results_rareprop$Method[overall_results_rareprop$Method == "RICE-CV"] <- "RICE-CV" 
overall_results_rareprop$Method[overall_results_rareprop$Method == "RICE-RV"] <- "RICE-RV" 
overall_results_rareprop$Method[overall_results_rareprop$Method == "LDPred"] <- "LDpred2"
overall_results_rareprop$Method[overall_results_rareprop$Method == "LASSOSum"] <- "Lassosum2"

overall_results_rareprop$Method1 <- overall_results_rareprop$Method
overall_results_rareprop$Method <- factor(overall_results_rareprop$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
overall_results_rareprop$Method1[overall_results_rareprop$Method1 == "RICE-RV"] <- "RICE-CV"
overall_results_rareprop$Method1 <- factor(overall_results_rareprop$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))
overall_results_rareprop$Ancestry <- factor(overall_results_rareprop$Ancestry,levels = c("AFR","AMR","EUR","SAS"))

results_CI$Method <- "RICE-RV"
overall_results <- left_join(overall_results,results_CI)
overall_results$Method <- factor(overall_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

overall_results$group1 <- "RICE-CV"
overall_results$group2 <- "RICE-CV"
overall_results$p.signif_beta <- ""
overall_results$p.signif_beta[overall_results$Method == "RICE-CV"] <- ifelse(overall_results$Beta_Adjusted_Lower_99[overall_results$Method == "RICE-RV"] > 0,"***",ifelse(overall_results$Beta_Adjusted_Lower_95[overall_results$Method == "RICE-RV"] > 0,"**",""))

overall_results$position[overall_results$Method == "RICE-CV"] <- overall_results$Beta[overall_results$Method == "RICE-CV"] + overall_results$Beta[overall_results$Method == "RICE-RV"] + 0.03
ylim_norm <- max(c(overall_results$Beta[overall_results$Method == "RICE-CV"] + overall_results$Beta[overall_results$Method == "RICE-RV"],overall_results$Beta)) + 0.05

results_rareprop_CI$Method <- "RICE-RV"
overall_results_rareprop <- left_join(overall_results_rareprop,results_rareprop_CI)
overall_results_rareprop$Method <- factor(overall_results_rareprop$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))

overall_results_rareprop$group1 <- "RICE-CV"
overall_results_rareprop$group2 <- "RICE-CV"
overall_results_rareprop$p.signif_beta <- ""
overall_results_rareprop$p.signif_beta[overall_results_rareprop$Method == "RICE-CV"] <- ifelse(overall_results_rareprop$Beta_Adjusted_Lower_99[overall_results_rareprop$Method == "RICE-RV"] > 0,"***",ifelse(overall_results_rareprop$Beta_Adjusted_Lower_95[overall_results_rareprop$Method == "RICE-RV"] > 0,"**",""))

overall_results_rareprop$position[overall_results_rareprop$Method == "RICE-CV"] <- overall_results_rareprop$Beta[overall_results_rareprop$Method == "RICE-CV"] + overall_results_rareprop$Beta[overall_results_rareprop$Method == "RICE-RV"] + 0.03
ylim_rareprop <- max(c(overall_results_rareprop$Beta[overall_results_rareprop$Method == "RICE-CV"] + overall_results_rareprop$Beta[overall_results_rareprop$Method == "RICE-RV"],overall_results_rareprop$Beta)) + 0.05

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

g1 <- ggplot(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_norm)+
  scale_fill_Publication()


ggsave(paste0("UKB_Simulation_Scaled_49173_Adjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)

g2 <- ggplot(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_norm)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Scaled_98343_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)

g3 <- ggplot(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_norm)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Unscaled_49173_Adjusted_Beta.png"),g3,width=10, height=6.18047,dpi = 300)

g4 <- ggplot(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_norm)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Unscaled_98343_Adjusted_Beta.png"),g4,width=10, height=6.18047,dpi = 300)

g1 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Scaled") & (overall_results_rareprop$Train_Size == "n = 49,173") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results_rareprop[(overall_results_rareprop$Scale == "Scaled") & (overall_results_rareprop$Train_Size == "n = 49,173") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_rareprop)+
  scale_fill_Publication()
ggsave(paste0("UKB_Simulation_RareProp_Scaled_49173_Adjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)

g2 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Scaled") & (overall_results_rareprop$Train_Size == "n = 98,343") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results_rareprop[(overall_results_rareprop$Scale == "Scaled") & (overall_results_rareprop$Train_Size == "n = 98,343") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_rareprop)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Scaled_98343_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)

g3 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Unscaled") & (overall_results_rareprop$Train_Size == "n = 49,173") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results_rareprop[(overall_results_rareprop$Scale == "Unscaled") & (overall_results_rareprop$Train_Size == "n = 49,173") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_rareprop) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Unscaled_49173_Adjusted_Beta.png"),g3,width=10, height=6.18047,dpi = 300)

g4 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Unscaled") & (overall_results_rareprop$Train_Size == "n = 98,343") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  stat_pvalue_manual(overall_results_rareprop[(overall_results_rareprop$Scale == "Unscaled") & (overall_results_rareprop$Train_Size == "n = 98,343") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),],
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2.5) +
  ylim(0,ylim_rareprop) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Unscaled_98343_Adjusted_Beta.png"),g4,width=10, height=6.18047,dpi = 300)

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#973999","#009cdf")), ...)
}

overall_results$Method <- as.character(overall_results$Method)
overall_results <- overall_results[overall_results$Method != "RICE-RV",]
overall_results$Method[overall_results$Method == "RICE-CV"] <- "RICE"
overall_results$Method <- factor(overall_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

overall_results_rareprop$Method <- as.character(overall_results_rareprop$Method)
overall_results_rareprop <- overall_results_rareprop[overall_results_rareprop$Method != "RICE-RV",]
overall_results_rareprop$Method[overall_results_rareprop$Method == "RICE-CV"] <- "RICE"
overall_results_rareprop$Method <- factor(overall_results_rareprop$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

overall_results <- overall_results[,c("Method","Scale","Causal_Prop","Train_Size","Ancestry","Beta","SE_Beta","R2","SE_R2")]

results_CI$Method <- "RICE"
overall_results <- left_join(overall_results,results_CI)
overall_results$Method <- factor(overall_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE"))

overall_results$group1 <- "RICE"
overall_results$group2 <- "RICE"
overall_results$p.signif_beta <- ""
overall_results$p.signif_beta[overall_results$Method == "RICE"] <- ifelse(overall_results$R2_Adjusted_Lower_99[overall_results$Method == "RICE"] > 0,"***",ifelse(overall_results$R2_Adjusted_Lower_95[overall_results$Method == "RICE"] > 0,"**",""))


g1 <- ggplot(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Scaled_49173_Adjusted_R2.png"),g1,width=10, height=6.18047,dpi = 300)

g2 <- ggplot(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Scaled_98343_Adjusted_R2.png"),g2,width=10, height=6.18047,dpi = 300)

g3 <- ggplot(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Unscaled_49173_Adjusted_R2.png"),g3,width=10, height=6.18047,dpi = 300)

g4 <- ggplot(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Unscaled_98343_Adjusted_R2.png"),g4,width=10, height=6.18047,dpi = 300)

g1 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Scaled") & (overall_results_rareprop$Train_Size == "n = 49,173") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Scaled_49173_Adjusted_R2.png"),g1,width=10, height=6.18047,dpi = 300)

g2 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Scaled") & (overall_results_rareprop$Train_Size == "n = 98,343") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Scaled_98343_Adjusted_R2.png"),g2,width=10, height=6.18047,dpi = 300)

g3 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Unscaled") & (overall_results_rareprop$Train_Size == "n = 49,173") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Unscaled_49173_Adjusted_R2.png"),g3,width=10, height=6.18047,dpi = 300)

g4 <- ggplot(overall_results_rareprop[(overall_results_rareprop$Scale == "Unscaled") & (overall_results_rareprop$Train_Size == "n = 98,343") & (overall_results_rareprop$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method, y=abs(R2),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("R2") + 
  theme_Publication() + 
  ylim(0,0.15) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_RareProp_Unscaled_98343_Adjusted_R2.png"),g4,width=10, height=6.18047,dpi = 300)
