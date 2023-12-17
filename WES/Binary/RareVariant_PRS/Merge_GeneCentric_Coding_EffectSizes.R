rm(list = ls())

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  for(i in 1:200){
    if(i == 1){
      coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_coding_sig_chr",i,".csv"))
    }else{
      coding_sig <- rbind(coding_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_coding_sig_chr",i,".csv")))
    }
  }
  coding_sig <- data.frame(coding_sig)
  coding_sig <- unique(coding_sig)
  write.csv(coding_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_Train_Effect_Sizes_All.csv"))
}