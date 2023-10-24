rm(list = ls())
for(i in 1:22){
  if(i == 1){
    sliding_window_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/sliding_window_sig_chr",i,".csv"))
  }else{
    sliding_window_sig <- rbind(sliding_window_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/sliding_window_sig_chr",i,".csv")))
  }
}
sliding_window_sig <- data.frame(sliding_window_sig)
sliding_window_sig <- unique(sliding_window_sig)
write.csv(sliding_window_sig,row.names = FALSE,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/SlidingWindow/Train_Effect_Sizes_All.csv")
