rm(list = ls())
for(i in 1:22){
  if(i == 1){
    coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/coding_sig_chr",i,".csv"))
  }else{
    coding_sig <- rbind(coding_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/coding_sig_chr",i,".csv")))
  }
}
coding_sig <- data.frame(coding_sig)
coding_sig <- unique(coding_sig)
write.csv(coding_sig,row.names = FALSE,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricCoding/Train_Effect_Sizes_All.csv")
