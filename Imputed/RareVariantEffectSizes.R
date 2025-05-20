rm(list = ls())
library("xlsx")
library(dplyr)

count <- 1
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  final_coef <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_Coefficients.csv"))
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  final_coef <- left_join(final_coef,coding_sig)
  final_coef <- final_coef[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef$Category[final_coef$Category == "synonymous"] <- "Synonymous"
  final_coef$Category[final_coef$Category == "plof"] <- "PLOF"
  final_coef$Category[final_coef$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef$Category[final_coef$Category == "missense"] <- "Missense"
  final_coef$Category[final_coef$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  final_coef$BETA[is.na(final_coef$BETA)] <- 0
  
  colnames(final_coef) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  if(count > 1){
    write.xlsx(final_coef, file = "UKB_Imputed.xlsx",sheetName = paste0(trait," UKB Imputed"), append = TRUE,row.names = FALSE)
  } else{
    write.xlsx(final_coef, file = "UKB_Imputed.xlsx",sheetName = paste0(trait," UKB Imputed"), append = FALSE)
  }
  count <- count + 1
}

rm(list = ls())
library("xlsx")
for(trait in c("Asthma","CAD","Breast","Prostate","T2D")){
  final_coef <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_Coefficients.csv"))
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  final_coef <- left_join(final_coef,coding_sig)
  final_coef <- final_coef[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef$Category[final_coef$Category == "synonymous"] <- "Synonymous"
  final_coef$Category[final_coef$Category == "plof"] <- "PLOF"
  final_coef$Category[final_coef$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef$Category[final_coef$Category == "missense"] <- "Missense"
  final_coef$Category[final_coef$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  final_coef$BETA[is.na(final_coef$BETA)] <- 0
  
  colnames(final_coef) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  write.xlsx(final_coef, file = "UKB_Imputed.xlsx",sheetName = paste0(trait," UKB Imputed"), append = TRUE)
}