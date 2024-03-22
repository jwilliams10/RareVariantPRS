rm(list = ls())

library(dplyr)

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  for(i in 1:200){
    if(i == 1){
      noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_noncoding_sigchr",i,".csv"))
    }else{
      noncoding_sig <- rbind(noncoding_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_noncoding_sigchr",i,".csv")))
    }
  }
  
  pvals <- unique(read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_noncoding_sig.csv")))
  
  noncoding_sig <- inner_join(noncoding_sig[,c("Gene","Chr","Category","Burden_Est")],pvals)
  noncoding_sig <- data.frame(noncoding_sig)
  noncoding_sig <- unique(noncoding_sig)
  write.csv(noncoding_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_Train_Effect_Sizes_All.csv"))
}

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
  
  pvals <- unique(read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_coding_sig.csv")))
  
  coding_sig <- inner_join(coding_sig[,c("Gene","Chr","Category","Burden_Est")],pvals)
  coding_sig <- data.frame(coding_sig)
  coding_sig <- unique(coding_sig)
  write.csv(coding_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_Train_Effect_Sizes_All.csv"))
}