rm(list = ls())

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  for(i in 1:200){
    if(i == 1){
      sliding_window_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_sliding_window_sig_array",i,".csv"))
    }else{
      sliding_window_sig <- rbind(sliding_window_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_sliding_window_sig_array",i,".csv")))
    }
  }
  sliding_window_sig <- data.frame(sliding_window_sig)
  sliding_window_sig <- unique(sliding_window_sig)
  write.csv(sliding_window_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_Train_Effect_Sizes_All.csv"))
}