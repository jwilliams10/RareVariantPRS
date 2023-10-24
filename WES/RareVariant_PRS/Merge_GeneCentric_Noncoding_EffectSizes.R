rm(list = ls())
for(i in 1:22){
  if(i == 1){
    noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/noncoding_sig_chr",i,".csv"))
  }else{
    noncoding_sig <- rbind(noncoding_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/noncoding_sig_chr",i,".csv")))
  }
}
noncoding_sig <- data.frame(noncoding_sig)
noncoding_sig <- unique(noncoding_sig)
write.csv(noncoding_sig,row.names = FALSE,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/GeneCentricNonCoding/Train_Effect_Sizes_All.csv")
