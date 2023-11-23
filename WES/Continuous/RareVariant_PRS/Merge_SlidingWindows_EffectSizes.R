rm(list = ls())

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  for(i in 1:22){
    if(i == 1){
      sliding_window_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_sliding_window_sig_chr",i,".csv"))
    }else{
      sliding_window_sig <- rbind(sliding_window_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_sliding_window_sig_chr",i,".csv")))
    }
  }
  sliding_window_sig <- data.frame(sliding_window_sig)
  sliding_window_sig <- unique(sliding_window_sig)
  write.csv(sliding_window_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/SlidingWindow/",trait,"_Train_Effect_Sizes_All.csv"))
}