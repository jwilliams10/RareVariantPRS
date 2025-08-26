rm(list = ls())

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

results_70 <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV_Alone <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Betas",i,".csv")))
  Best_Betas_RICECV_Alone$Method <- "RICE-CV-Alone"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  Best_Betas_RICECV_Alone <- Best_Betas_RICECV_Alone[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV_Alone,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_70 <- rbind(results_70,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","i","Y_train","index_mat")))
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







load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_35 <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV_Alone <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Betas",i,".csv")))
  Best_Betas_RICECV_Alone$Method <- "RICE-CV-Alone"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  Best_Betas_RICECV_Alone <- Best_Betas_RICECV_Alone[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV_Alone,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_35 <- rbind(results_35,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_35","i","Y_train","index_mat")))
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





load("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_rareprop_70 <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV_Alone <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Combined_Common_PRS/Best_Betas",i,".csv")))
  Best_Betas_RICECV_Alone$Method <- "RICE-CV-Alone"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  Best_Betas_RICECV_Alone <- Best_Betas_RICECV_Alone[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV_Alone,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_rareprop_70 <- rbind(results_rareprop_70,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_35","results_rareprop_70","i","Y_train","index_mat")))
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





load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_rareprop_35 <- NULL

for(i in 1:length(Y_train)){
  
  Best_Betas_CT <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/Best_Betas",i,".csv")))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/LDPred2/Best_Betas",i,".csv")))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/LASSOSUM2/Best_Betas",i,".csv")))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV_Alone <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/Combined_Common_PRS/Best_Betas",i,".csv")))
  Best_Betas_RICECV_Alone$Method <- "RICE-CV-Alone"
  
  Best_Betas_RICECV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv")))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv")))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  Best_Betas_CT <- Best_Betas_CT[,colnames(Best_Betas_RICECV)]
  Best_Betas_LDPred <- Best_Betas_LDPred[,colnames(Best_Betas_RICECV)]
  Best_Betas_LASSOSum <- Best_Betas_LASSOSum[,colnames(Best_Betas_RICECV)]
  Best_Betas_RICECV_Alone <- Best_Betas_RICECV_Alone[,colnames(Best_Betas_RICECV)]
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV_Alone,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_rareprop_35 <- rbind(results_rareprop_35,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_70","results_35","results_rareprop_70","results_rareprop_35","i","Y_train","index_mat")))
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






results <- rbind(results_35,results_70)
results_rareprop <- rbind(results_rareprop_35,results_rareprop_70)

results$Train_Size <- format(results$Train_Size,big.mark=",", trim=TRUE)
results_rareprop$Train_Size <- format(results_rareprop$Train_Size,big.mark=",", trim=TRUE)

results$Train_Size <- paste0("n = ",results$Train_Size)
results_rareprop$Train_Size <- paste0("n = ",results_rareprop$Train_Size)

rm(list=setdiff(ls(), c("results","results_rareprop")))

results$Beta[results$Method %in% c("LDPred","LASSOSum")] <- -1*results$Beta[results$Method %in% c("LDPred","LASSOSum")]
results$Beta[results$Beta < 0] <- 0
results$Beta[results$Beta > 1] <- 0

results_rareprop$Beta[results_rareprop$Method %in% c("LDPred","LASSOSum")] <- -1*results_rareprop$Beta[results_rareprop$Method %in% c("LDPred","LASSOSum")]
results_rareprop$Beta[results_rareprop$Beta < 0] <- 0
results_rareprop$Beta[results_rareprop$Beta > 1] <- 0

results <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results,mean)
results_rareprop <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results_rareprop,mean)

overall_results <- results[results$Method %in% c("CT","LASSOSum","LDPred","RICE-CV-Alone","RICE-CV","RICE-RV"),]
overall_results_rareprop <- results_rareprop[results_rareprop$Method %in% c("CT","LASSOSum","LDPred","RICE-CV-Alone","RICE-CV","RICE-RV"),]

overall_results <- overall_results[overall_results$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results$Method[overall_results$Method == "RICE-CV"] <- "RICE-CV" 
overall_results$Method[overall_results$Method == "RICE-RV"] <- "RICE-RV" 
overall_results$Method[overall_results$Method == "LDPred"] <- "LDpred2"
overall_results$Method[overall_results$Method == "LASSOSum"] <- "Lassosum2"

overall_results$Method1 <- overall_results$Method
overall_results$Method <- factor(overall_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV-Alone","RICE-CV"))
overall_results$Ancestry <- factor(overall_results$Ancestry,levels = c("AFR","AMR","EUR","SAS"))

overall_results_rareprop <- overall_results_rareprop[overall_results_rareprop$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results_rareprop$Method[overall_results_rareprop$Method == "RICE-CV"] <- "RICE-CV" 
overall_results_rareprop$Method[overall_results_rareprop$Method == "RICE-RV"] <- "RICE-RV" 
overall_results_rareprop$Method[overall_results_rareprop$Method == "LDPred"] <- "LDpred2"
overall_results_rareprop$Method[overall_results_rareprop$Method == "LASSOSum"] <- "Lassosum2"

overall_results_rareprop$Method1 <- overall_results_rareprop$Method
overall_results_rareprop$Method <- factor(overall_results_rareprop$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV-Alone","RICE-CV"))
overall_results_rareprop$Ancestry <- factor(overall_results_rareprop$Ancestry,levels = c("AFR","AMR","EUR","SAS"))

overall_results$Causality_Structure <- "All"
overall_results_rareprop$Causality_Structure <- "Proportion"
overall_results <- rbind(overall_results,overall_results_rareprop)

beta_ratios <- NULL
for(caus_struct in c("All","Proportion")){
  for(i in c("Scaled","Unscaled")){
    for(k in c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")){
      for(q in c("n = 49,173","n = 98,343")){
        for(ancs in c("AFR","AMR","EUR","SAS")){
          data <- overall_results[overall_results$Scale == i & overall_results$Causal_Prop == k & overall_results$Train_Size == q & overall_results$Ancestry == ancs & overall_results$Causality_Structure == caus_struct,]
          
          ratio_RV_Add_R2 <- (data$R2[data$Method == "RICE-CV"])/data$R2[data$Method == "RICE-CV-Alone"] - 1
          ratio_CV_Alone_Best_R2 <- data$R2[data$Method == "RICE-CV-Alone"]/max(data$R2[!(data$Method %in% c("RICE-CV-Alone","RICE-CV","RICE-RV"))]) - 1
          
          beta_ratios <- rbind(beta_ratios,data.frame(Causality_Structure = caus_struct,Scale = i, Causal_Prop = k, Training_Size = q, Ancestry = ancs,ratio_RV_Add_R2 =ratio_RV_Add_R2,ratio_CV_Alone_Best_R2 = ratio_CV_Alone_Best_R2))
          
        }
      }
    }
  } 
}

mean(beta_ratios$ratio_RV_Add_R2[beta_ratios$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")])

mean(beta_ratios$ratio_CV_Alone_Best_R2[beta_ratios$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")])

