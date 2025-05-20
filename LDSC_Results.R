rm(list = ls())
library(readr)
library(stringr)

LDSC_Results <- NULL

for(trait in c("Asthma","CAD","T2D","Breast","Prostate","BMI","Height","HDL","LDL","logTG","TC")){

  if(trait %in% c("BMI","Height","HDL","LDL","logTG","TC")){
    LDSC_WES_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Results.log"), skip = 20)
    LDSC_WES_Results <- as.data.frame(LDSC_WES_Results)
    LDSC_WES_Intercept <- gsub("Intercept: ","",LDSC_WES_Results[str_detect(LDSC_WES_Results[,1],"Intercept:"),1])
    LDSC_WES_h2 <- gsub("Total Observed scale h2: ","",LDSC_WES_Results[str_detect(LDSC_WES_Results[,1],"Total Observed scale h2:"),1]) 
  }else{
    LDSC_WES_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/LDSC/",trait,"_LDSC_Results.log"), skip = 20)
    LDSC_WES_Results <- as.data.frame(LDSC_WES_Results)
    LDSC_WES_Intercept <- gsub("Intercept: ","",LDSC_WES_Results[str_detect(LDSC_WES_Results[,1],"Intercept:"),1])
    LDSC_WES_h2 <- gsub("Total Observed scale h2: ","",LDSC_WES_Results[str_detect(LDSC_WES_Results[,1],"Total Observed scale h2:"),1]) 
  }
  
  LDSC_Imputed_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDSC/",trait,"_LDSC_Imputed_Results.log"), skip = 20)
  LDSC_Imputed_Results <- as.data.frame(LDSC_Imputed_Results)
  LDSC_Imputed_Intercept <- gsub("Intercept: ","",LDSC_Imputed_Results[str_detect(LDSC_Imputed_Results[,1],"Intercept:"),1])
  LDSC_Imputed_h2 <- gsub("Total Observed scale h2: ","",LDSC_Imputed_Results[str_detect(LDSC_Imputed_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_WGS_Results <- read_csv(paste0("/data/williamsjacr/Clean_UKB_WGS_Sumstats/LDSC/",trait,"_LDSC_Results.log"), skip = 20)
  LDSC_WGS_Results <- as.data.frame(LDSC_WGS_Results)
  LDSC_WGS_Intercept <- gsub("Intercept: ","",LDSC_WGS_Results[str_detect(LDSC_WGS_Results[,1],"Intercept:"),1])
  LDSC_WGS_h2 <- gsub("Total Observed scale h2: ","",LDSC_WGS_Results[str_detect(LDSC_WGS_Results[,1],"Total Observed scale h2:"),1])
  
  if(trait %in% c("BMI","Height","HDL","LDL","logTG","TC")){
    LDSC_EUR_AoU_Results <- read_csv(paste0("/data/williamsjacr/Cleaned_AoU_SumStats/LDSC/EUR_",trait,"_LDSC_Results.log"), skip = 20)
    LDSC_EUR_AoU_Results <- as.data.frame(LDSC_EUR_AoU_Results)
    LDSC_EUR_AoU_Intercept <- gsub("Intercept: ","",LDSC_EUR_AoU_Results[str_detect(LDSC_EUR_AoU_Results[,1],"Intercept:"),1])
    LDSC_EUR_AoU_h2 <- gsub("Total Observed scale h2: ","",LDSC_EUR_AoU_Results[str_detect(LDSC_EUR_AoU_Results[,1],"Total Observed scale h2:"),1])
    
    LDSC_AFR_AoU_Results <- read_csv(paste0("/data/williamsjacr/Cleaned_AoU_SumStats/LDSC/AFR_",trait,"_LDSC_Results.log"), skip = 20)
    LDSC_AFR_AoU_Results <- as.data.frame(LDSC_AFR_AoU_Results)
    LDSC_AFR_AoU_Intercept <- gsub("Intercept: ","",LDSC_AFR_AoU_Results[str_detect(LDSC_AFR_AoU_Results[,1],"Intercept:"),1])
    LDSC_AFR_AoU_h2 <- gsub("Total Observed scale h2: ","",LDSC_AFR_AoU_Results[str_detect(LDSC_AFR_AoU_Results[,1],"Total Observed scale h2:"),1])
    
    LDSC_AMR_AoU_Results <- read_csv(paste0("/data/williamsjacr/Cleaned_AoU_SumStats/LDSC/AMR_",trait,"_LDSC_Results.log"), skip = 20)
    LDSC_AMR_AoU_Results <- as.data.frame(LDSC_AMR_AoU_Results)
    LDSC_AMR_AoU_Intercept <- gsub("Intercept: ","",LDSC_AMR_AoU_Results[str_detect(LDSC_AMR_AoU_Results[,1],"Intercept:"),1])
    LDSC_AMR_AoU_h2 <- gsub("Total Observed scale h2: ","",LDSC_AMR_AoU_Results[str_detect(LDSC_AMR_AoU_Results[,1],"Total Observed scale h2:"),1]) 
    
    LDSC_Results <- rbind(LDSC_Results,data.frame(trait = trait,Ancestry = c("EUR","EUR","EUR","AFR","AMR","EUR"),Source = c("UKB WES","UKB Imputed + WES","UKB WGS","AoU","AoU","AoU"),
                                                  h2 = c(LDSC_WES_h2,LDSC_Imputed_h2,LDSC_WGS_h2,LDSC_AFR_AoU_h2,LDSC_AMR_AoU_h2,LDSC_EUR_AoU_h2),
                                                  Intercept = c(LDSC_WES_Intercept,LDSC_Imputed_Intercept,LDSC_WGS_Intercept,LDSC_AFR_AoU_Intercept,LDSC_AMR_AoU_Intercept,LDSC_EUR_AoU_Intercept)))
  }else{
    LDSC_Results <- rbind(LDSC_Results,data.frame(trait = trait,Ancestry = c("EUR","EUR","EUR"),Source = c("UKB WES","UKB Imputed + WES","UKB WGS"),
                                                  h2 = c(LDSC_WES_h2,LDSC_Imputed_h2,LDSC_WGS_h2),
                                                  Intercept = c(LDSC_WES_Intercept,LDSC_Imputed_Intercept,LDSC_WGS_Intercept)))
  }
  
}

write.csv(LDSC_Results,file = "LDSC_Results.csv",row.names = FALSE)