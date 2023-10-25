rm(list = ls())

library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")

## Match Columns 

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/SlidingWindow/chr22.Rdata")
slidingwindow_sig <- data.frame(Chr = unlist(results_sliding_window[,1]),Start = unlist(results_sliding_window[,2]),End = unlist(results_sliding_window[,3]))
slidingwindow_sig$Chr <- as.numeric(slidingwindow_sig$Chr)
slidingwindow_sig_overall <- slidingwindow_sig

slidingwindow_sig <- read.csv("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes1.csv")
slidingwindow_sig <- left_join(slidingwindow_sig_overall,slidingwindow_sig)

## Match Rows

obj_nullmodel_full <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

obj_nullmodel_sub <- get(load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Validation_Null_Model1.RData"))

## Filter

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/SlidingWindow.RData")

G_star_sliding_window <- G_star_sliding_window[which(obj_nullmodel_full$id_include %in% obj_nullmodel_sub$id_include),]
G_star_sliding_window <- G_star_sliding_window[,!is.na(slidingwindow_sig$Burden_1_1)]


### PRS

thresholds <- c(1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,5e-04,1e-03,5e-03,1e-02,5e-02,1e-01,5e-01,1.0)

for(i in 1:length(Y_validation)){
  slidingwindow_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/Train_EffectSizes",i,".csv"))
  prs_mat <- NULL
  
  for(j in 1:2){
    if(j == 1){
      p_values <- slidingwindow_sig$Burden_pvalue
    }else{
      p_values <- slidingwindow_sig$STAAR_O
    }
    for(k in 1:length(thresholds)){
      beta <- matrix(0,ncol = 1,nrow = ncol(G_star_sliding_window))
      beta[p_values < thresholds[k]] <- slidingwindow_sig$Burden_Est[p_values < thresholds[k]]
      prs_mat <- cbind(prs_mat,G_star_sliding_window%*%beta)
    }
  }
  save(prs_mat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/validation_prs_mat",i,".RData"))
}


