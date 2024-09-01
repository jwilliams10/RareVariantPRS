
rm(list = ls())
library("xlsx")
count <- 1
for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  final_coef <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"))
  final_coef$Category[final_coef$Category == "synonymous"] <- "Synonymous"
  final_coef$Category[final_coef$Category == "plof"] <- "PLOF"
  final_coef$Category[final_coef$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef$Category[final_coef$Category == "missense"] <- "Missense"
  final_coef$Category[final_coef$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  colnames(final_coef) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  if(count > 1){
    write.xlsx(final_coef, file = "UKB_WES.xlsx",sheetName = paste0(trait," UKB WES"), append = TRUE,row.names = FALSE)
  } else{
    write.xlsx(final_coef, file = "UKB_WES.xlsx",sheetName = paste0(trait," UKB WES"), append = FALSE)
  }
  count <- count + 1
}

rm(list = ls())
library("xlsx")
for(trait in c("Asthma","CAD","Breast","Prostate","T2D")){
  final_coef <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"))
  final_coef$Category[final_coef$Category == "synonymous"] <- "Synonymous"
  final_coef$Category[final_coef$Category == "plof"] <- "PLOF"
  final_coef$Category[final_coef$Category == "plof_ds"] <- "PLOF (DS)"
  final_coef$Category[final_coef$Category == "missense"] <- "Missense"
  final_coef$Category[final_coef$Category == "disruptive_missense"] <- "Disruptive Missense"
  
  colnames(final_coef) <- c("Gene","Chromosome","Category","Number of SNV","Burden(1,1)","STAAR-B(1,1)","Estimated Burden Effect Size")
  write.xlsx(final_coef, file = "UKB_WES.xlsx",sheetName = paste0(trait," UKB WES"), append = TRUE)
}