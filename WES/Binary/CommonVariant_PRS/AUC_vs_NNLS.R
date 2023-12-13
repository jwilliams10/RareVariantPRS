rm(list = ls())

a <- data.frame(trait = NULL, anc = NULL, diff = NULL)

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  for(anc in c("Eur","NonEur","UNK","MIX","AFR","SAS","EAS")){
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_",anc,"_NNLS.RData"))
    if(sum(is.na(SL.result))>0){
      b1 <- NA
    }else{
      b1 <- SL.result$AUC  
    }
    
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_sl_result_All_",anc,".RData"))
    if(sum(is.na(SL.result))>0){
      b2 <- NA
    }else{
      b2 <- SL.result$AUC
    }
    
    a <- rbind(a, data.frame(trait = trait, anc = anc, diff = b2 - b1))
  }
}


load("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/Asthma_sl_result_All_AFR_NNLS.RData")
SL.result$AUC
load("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/Asthma_sl_result_All_AFR.RData")
SL.result$AUC
