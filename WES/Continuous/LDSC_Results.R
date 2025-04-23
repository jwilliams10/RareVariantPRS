rm(list = ls())
library(readr)
library(stringr)

LDSC_Results <- NULL

for(trait in c("BMI","Height","HDL","LDL","logTG","TC")){
  LDSC_20PCs_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_20PCs_LDSC_Results.log"), skip = 20)
  LDSC_20PCs_Results <- as.data.frame(LDSC_20PCs_Results)
  LDSC_20PCs_Intercept <- gsub("Intercept: ","",LDSC_20PCs_Results[str_detect(LDSC_20PCs_Results[,1],"Intercept:"),1])
  LDSC_20PCs_h2 <- gsub("Total Observed scale h2: ","",LDSC_20PCs_Results[str_detect(LDSC_20PCs_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_40PCs_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_40PCs_LDSC_Results.log"), skip = 20)
  LDSC_40PCs_Results <- as.data.frame(LDSC_40PCs_Results)
  LDSC_40PCs_Intercept <- gsub("Intercept: ","",LDSC_40PCs_Results[str_detect(LDSC_40PCs_Results[,1],"Intercept:"),1])
  LDSC_40PCs_h2 <- gsub("Total Observed scale h2: ","",LDSC_40PCs_Results[str_detect(LDSC_40PCs_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Original_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_OriginalPlink_LDSC_Results.log"), skip = 20)
  LDSC_Original_Results <- as.data.frame(LDSC_Original_Results)
  LDSC_Original_Intercept <- gsub("Intercept: ","",LDSC_Original_Results[str_detect(LDSC_Original_Results[,1],"Intercept:"),1])
  LDSC_Original_h2 <- gsub("Total Observed scale h2: ","",LDSC_Original_Results[str_detect(LDSC_Original_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Regenie_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Results.log"), skip = 20)
  LDSC_Regenie_Results <- as.data.frame(LDSC_Regenie_Results)
  LDSC_Regenie_Intercept <- gsub("Intercept: ","",LDSC_Regenie_Results[str_detect(LDSC_Regenie_Results[,1],"Intercept:"),1])
  LDSC_Regenie_h2 <- gsub("Total Observed scale h2: ","",LDSC_Regenie_Results[str_detect(LDSC_Regenie_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_RankNormal_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"adj_norm_ranknormal_LDSC_Results.log"), skip = 20)
  LDSC_RankNormal_Results <- as.data.frame(LDSC_RankNormal_Results)
  LDSC_RankNormal_Intercept <- gsub("Intercept: ","",LDSC_RankNormal_Results[str_detect(LDSC_RankNormal_Results[,1],"Intercept:"),1])
  LDSC_RankNormal_h2 <- gsub("Total Observed scale h2: ","",LDSC_RankNormal_Results[str_detect(LDSC_RankNormal_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Results <- rbind(LDSC_Results,data.frame(trait = trait,method = c("regenie 10 PCs","plink 10 PCs","plink 20 PCs","plink 40 PCs","plink rank normal 10 PCs"),
                                                h2 = c(LDSC_Regenie_h2,LDSC_Original_h2,LDSC_20PCs_h2,LDSC_40PCs_h2,LDSC_RankNormal_h2),
                                                Intercept = c(LDSC_Regenie_Intercept,LDSC_Original_Intercept,LDSC_20PCs_Intercept,LDSC_40PCs_Intercept,LDSC_RankNormal_Intercept)))
  
}

write.csv(LDSC_Results,file = "LDSC_Results.csv",row.names = FALSE)