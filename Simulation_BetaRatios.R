rm(list = ls())

library(dplyr)

index_mat <- NULL

causalprop_vec <- c(0.2,0.05,0.01,0.001,0.0005)
scale <- c(0,1)

count <- 1

for(j in 1:length(causalprop_vec)){
  for(q in 1:length(scale)){
    for(l in 1:20){
      index_mat <- rbind(index_mat,data.frame(i = count,Causal_Prop = causalprop_vec[j],Scale = scale[q]))
      count <- count + 1
    }
  }
}

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_70 <- NULL

for(i in 1:length(Y_train)){
  Best_Betas_CV_SL <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/Best_Betas",i,".csv"))
  Best_Betas_CV_SL$Method <- "CV_SL"
  
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RV_SL <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_Betas",i,".csv"))
  Best_Betas_RV_SL$Method <- "RV_SL"
  
  Best_Betas_lm <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Best_Betas",i,".csv"))
  
  if(sum(Best_Betas_CV_SL$beta_raw == Best_Betas_LDPred$beta_raw) == 7){
    Best_Betas_CV_SL$beta_raw <- -1*Best_Betas_CV_SL$beta_raw
    Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"] <- -1*Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"]
  }
  
  if(sum(Best_Betas_CV_SL$beta_raw == Best_Betas_LASSOSum$beta_raw) == 7){
    Best_Betas_CV_SL$beta_raw <- -1*Best_Betas_CV_SL$beta_raw
    Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"] <- -1*Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"]
  }
  
  betas_tmp <- rbind(Best_Betas_CV_SL,Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RV_SL,Best_Betas_lm)
  
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

results_70 <- data.frame(Scale = results_70$Scale, Causal_Prop = results_70$Causal_Prop, Method = results_70$Method,Ancestry = results_70$ancestry,Beta = results_70$beta_raw,SE_Beta = results_70$se_raw)
results_70$Train_Size <- nrow(Y_train[[1]])







load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_35 <- NULL

for(i in 1:length(Y_train)){
  Best_Betas_CV_SL <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/Best_Betas",i,".csv"))
  Best_Betas_CV_SL$Method <- "CV_SL"
  
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RV_SL <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/Best_Betas",i,".csv"))
  Best_Betas_RV_SL$Method <- "RV_SL"
  
  Best_Betas_lm <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Best_Betas",i,".csv"))
  
  if(sum(Best_Betas_CV_SL$beta_raw == Best_Betas_LDPred$beta_raw) == 7){
    Best_Betas_CV_SL$beta_raw <- -1*Best_Betas_CV_SL$beta_raw
    Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"] <- -1*Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"]
  }
  
  if(sum(Best_Betas_CV_SL$beta_raw == Best_Betas_LASSOSum$beta_raw) == 7){
    Best_Betas_CV_SL$beta_raw <- -1*Best_Betas_CV_SL$beta_raw
    Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"] <- -1*Best_Betas_lm$beta_raw[Best_Betas_lm$Method == "CV"]
  }
  
  betas_tmp <- rbind(Best_Betas_CV_SL,Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RV_SL,Best_Betas_lm)
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

results_35 <- data.frame(Scale = results_35$Scale, Causal_Prop = results_35$Causal_Prop, Method = results_35$Method,Ancestry = results_35$ancestry,Beta = results_35$beta_raw,SE_Beta = results_35$se_raw)
results_35$Train_Size <- nrow(Y_train[[1]])






results <- rbind(results_35,results_70)
results$Train_Size <- format(results$Train_Size,big.mark=",", trim=TRUE)

results$Train_Size <- paste0("n = ",results$Train_Size)

rm(list=setdiff(ls(), c("results")))

results$Beta[results$Beta < 0 & results$Method %in% c("LDPred","LASSOSum")] <- -1*results$Beta[results$Beta < 0 & results$Method %in% c("LDPred","LASSOSum")]
results$Beta[results$Beta < 0] <- 0

results <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results,mean)
overall_results <- results[results$Method %in% c("CT","LASSOSum","LDPred","CV_SL","CV","RV"),]

# overall_results <- read.csv("Desktop/RareVariantPRS_Results/Overall_Results_SimStudy.csv")

overall_results <- overall_results[overall_results$Ancestry %in% c("AFR","EUR","SAS","AMR","EAS"),]
overall_results <- overall_results[overall_results$Method %in% c("CT","LASSOSum","LDPred","CV","RV"),]
overall_results$Method[overall_results$Method == "CV"] <- "RICE-CV" 
overall_results$Method[overall_results$Method == "RV"] <- "RICE-RV" 
overall_results$Method[overall_results$Method == "LDPred"] <- "LDpred2"
overall_results$Method[overall_results$Method == "LASSOSum"] <- "Lassosum2"

overall_results <- overall_results[overall_results$Ancestry %in% c("AFR","AMR","EUR","SAS"),]


beta_ratios <- NULL
for(i in c("Scaled","Unscaled")){
  for(k in c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")){
    for(q in c("n = 49,173","n = 98,343")){
      for(ancs in c("AFR","AMR","EUR","SAS")){
        data <- overall_results[overall_results$Scale == i & overall_results$Causal_Prop == k & overall_results$Train_Size == q & overall_results$Ancestry == ancs,]
        
        ratio_CV_RV_BestofRest <- (data$Beta[data$Method == "RICE-CV"] + data$Beta[data$Method == "RICE-RV"])/max(data$Beta[!(data$Method %in% c("RICE-CV","RICE-RV"))]) - 1
        ratio_RV_CV <- data$Beta[data$Method == "RICE-RV"]/data$Beta[data$Method == "RICE-CV"]
        ratio_RV_Total <- data$Beta[data$Method == "RICE-RV"]/(data$Beta[data$Method == "RICE-CV"] + data$Beta[data$Method == "RICE-RV"])
        
        beta_ratios <- rbind(beta_ratios,data.frame(Scale = i, Causal_Prop = k, Training_Size = q, Ancestry = ancs,Ratio_RICECV_RV_Best = ratio_CV_RV_BestofRest,Ratio_RV_CV = ratio_RV_CV,Ratio_RV_Total = ratio_RV_Total))
      }
    }
  }
}

mean(beta_ratios$Ratio_RV_Total[beta_ratios$Ancestry == "AFR" & beta_ratios$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")])
mean(beta_ratios$Ratio_RV_Total[beta_ratios$Ancestry == "EUR" & beta_ratios$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")])
mean(beta_ratios$Ratio_RV_Total[beta_ratios$Ancestry == "SAS" & beta_ratios$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")])
mean(beta_ratios$Ratio_RV_Total[beta_ratios$Ancestry == "AMR" & beta_ratios$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")])


