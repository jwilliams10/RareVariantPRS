rm(list = ls())
library("xlsx")
library(dplyr)
count <- 1
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  final_coef <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"))
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
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
    write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB WES"), append = TRUE,row.names = FALSE)
  } else{
    write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB WES"), append = FALSE)
  }
  count <- count + 1
}

rm(list = ls())
library("xlsx")
for(trait in c("Asthma","CAD","Breast","Prostate","T2D")){
  final_coef <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"))
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  final_coef <- left_join(final_coef,coding_sig)
  final_coef <- final_coef[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef$Category[final_coef$Category == "synonymous"] <- "Synonymous"
  final_coef$Category[final_coef$Category == "plof"] <- "PLOF"
  final_coef$Category[final_coef$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef$Category[final_coef$Category == "missense"] <- "Missense"
  final_coef$Category[final_coef$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  final_coef$BETA[is.na(final_coef$BETA)] <- 0
  
  colnames(final_coef) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB WES"), append = TRUE,row.names = FALSE)
}

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
  write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB Imputed"), append = TRUE,row.names = FALSE)
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
  write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB Imputed"), append = TRUE,row.names = FALSE)
}

rm(list = ls())
library("xlsx")
library(dplyr)
count <- 1

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  final_coef_coding <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"_Coding_Final_Coefficients.csv"))
  coding_sig <- read.csv("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/coding_sig.csv")
  coding_sig <- coding_sig[coding_sig$Trait == trait,]
  coding_sig <- subset(coding_sig,select = -Trait)
  final_coef_coding <- left_join(final_coef_coding,coding_sig)
  final_coef_coding <- final_coef_coding[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef_coding$Category[final_coef_coding$Category == "synonymous"] <- "Synonymous"
  final_coef_coding$Category[final_coef_coding$Category == "plof"] <- "PLOF"
  final_coef_coding$Category[final_coef_coding$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef_coding$Category[final_coef_coding$Category == "missense"] <- "Missense"
  final_coef_coding$Category[final_coef_coding$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  final_coef_coding$BETA[is.na(final_coef_coding$BETA)] <- 0
  
  final_coef_noncoding <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"_Noncoding_Final_Coefficients.csv"))
  coding_sig <- read.csv("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/noncoding_sig.csv")
  coding_sig <- coding_sig[coding_sig$Trait == trait,]
  coding_sig <- subset(coding_sig,select = -Trait)
  final_coef_noncoding <- left_join(final_coef_noncoding,coding_sig)
  final_coef_noncoding <- final_coef_noncoding[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef_noncoding$Category[final_coef_noncoding$Category == "upstream"] <- "Upstream"  
  final_coef_noncoding$Category[final_coef_noncoding$Category == "downstream"] <- "Downstream"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "enhancer_CAGE"] <- "Enhancer (CAGE)"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "promoter_DHS"] <- "Promoter (DHS)"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "enhancer_DHS"] <- "Enhancer (DHS)"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "promoter_CAGE"] <- "Promoter (CAGE)"
  
  final_coef_noncoding$BETA[is.na(final_coef_noncoding$BETA)] <- 0
  
  colnames(final_coef_coding) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  colnames(final_coef_noncoding) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  final_coef <- rbind(final_coef_coding,final_coef_noncoding)

  write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB WGS"), append = TRUE,row.names = FALSE)
  count <- count + 1
}

rm(list = ls())
library("xlsx")
for(trait in c("Asthma","CAD","Breast","Prostate","T2D")){
  final_coef_coding <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/",trait,"_Coding_Final_Coefficients.csv"))
  coding_sig <- read.csv("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/coding_sig.csv")
  coding_sig <- coding_sig[coding_sig$Trait == trait,]
  coding_sig <- subset(coding_sig,select = -Trait)
  final_coef_coding <- left_join(final_coef_coding,coding_sig)
  final_coef_coding <- final_coef_coding[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef_coding$Category[final_coef_coding$Category == "synonymous"] <- "Synonymous"
  final_coef_coding$Category[final_coef_coding$Category == "plof"] <- "PLOF"
  final_coef_coding$Category[final_coef_coding$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef_coding$Category[final_coef_coding$Category == "missense"] <- "Missense"
  final_coef_coding$Category[final_coef_coding$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  final_coef_coding$BETA[is.na(final_coef_coding$BETA)] <- 0
  
  final_coef_noncoding <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/",trait,"_Noncoding_Final_Coefficients.csv"))
  coding_sig <- read.csv("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/noncoding_sig.csv")
  coding_sig <- coding_sig[coding_sig$Trait == trait,]
  coding_sig <- subset(coding_sig,select = -Trait)
  final_coef_noncoding <- left_join(final_coef_noncoding,coding_sig)
  final_coef_noncoding <- final_coef_noncoding[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","BETA")]
  
  final_coef_noncoding$Category[final_coef_noncoding$Category == "upstream"] <- "Upstream"  
  final_coef_noncoding$Category[final_coef_noncoding$Category == "downstream"] <- "Downstream"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "enhancer_CAGE"] <- "Enhancer (CAGE)"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "promoter_DHS"] <- "Promoter (DHS)"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "enhancer_DHS"] <- "Enhancer (DHS)"
  final_coef_noncoding$Category[final_coef_noncoding$Category == "promoter_CAGE"] <- "Promoter (CAGE)"
  
  final_coef_noncoding$BETA[is.na(final_coef_noncoding$BETA)] <- 0
  
  colnames(final_coef_coding) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  colnames(final_coef_noncoding) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  final_coef <- rbind(final_coef_coding,final_coef_noncoding)
  write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," UKB WGS"), append = TRUE,row.names = FALSE)
}


rm(list = ls())
library("xlsx")
library(dplyr)

count <- 1
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  final_coef <- read.csv(paste0("/data/williamsjacr/AoU_Results/",trait,"_final_coef.csv"))
  final_coef <- final_coef[,c("Gene","Chr","Category","Number_SNV","Burden_1_1","STAARB","Beta")]
  
  final_coef$Category[final_coef$Category == "synonymous"] <- "Synonymous"
  final_coef$Category[final_coef$Category == "plof"] <- "PLOF"
  final_coef$Category[final_coef$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef$Category[final_coef$Category == "missense"] <- "Missense"
  final_coef$Category[final_coef$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  final_coef$Beta[is.na(final_coef$Beta)] <- 0
  
  colnames(final_coef) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  write.xlsx(final_coef, file = "Supplementary_Data.xlsx",sheetName = paste0(trait," AoU"), append = TRUE,row.names = FALSE)

  count <- count + 1
}