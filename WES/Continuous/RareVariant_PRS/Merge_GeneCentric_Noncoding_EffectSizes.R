rm(list = ls())

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  for(i in 1:22){
    if(i == 1){
      noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_noncoding_sigchr",i,".csv"))
    }else{
      noncoding_sig <- rbind(noncoding_sig, read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_noncoding_sigchr",i,".csv")))
    }
  }
  noncoding_sig <- data.frame(noncoding_sig)
  noncoding_sig <- unique(noncoding_sig)
  write.csv(noncoding_sig,row.names = FALSE,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_Train_Effect_Sizes_All.csv"))
}