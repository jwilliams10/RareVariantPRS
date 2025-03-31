rm(list = ls())
library(readr)
library(stringr)

LDSC_Results <- NULL

for(trait in c("BMI","HDL","LDL","logTG","TC","Height")){
  LDSC_20PCs_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_20PCs_LDSC_Results.log"), skip = 20)
  LDSC_20PCs_Results <- as.data.frame(LDSC_20PCs_Results)
  LDSC_20PCs_Intercept <- gsub("Intercept: ","",LDSC_20PCs_Results[str_detect(LDSC_20PCs_Results[,1],"Intercept:"),1])
  LDSC_20PCs_h2 <- gsub("Total Observed scale h2: ","",LDSC_20PCs_Results[str_detect(LDSC_20PCs_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Original_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_OriginalPlink_LDSC_Results.log"), skip = 20)
  LDSC_Original_Results <- as.data.frame(LDSC_Original_Results)
  LDSC_Original_Intercept <- gsub("Intercept: ","",LDSC_Original_Results[str_detect(LDSC_Original_Results[,1],"Intercept:"),1])
  LDSC_Original_h2 <- gsub("Total Observed scale h2: ","",LDSC_Original_Results[str_detect(LDSC_Original_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Regenie_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Results.log"), skip = 20)
  LDSC_Regenie_Results <- as.data.frame(LDSC_Regenie_Results)
  LDSC_Regenie_Intercept <- gsub("Intercept: ","",LDSC_Regenie_Results[str_detect(LDSC_Regenie_Results[,1],"Intercept:"),1])
  LDSC_Regenie_h2 <- gsub("Total Observed scale h2: ","",LDSC_Regenie_Results[str_detect(LDSC_Regenie_Results[,1],"Total Observed scale h2:"),1])
  
  if(trait %in% c("LDL","TC")){
    LDSC_Consortium_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_statin_adj_LDSC_Consortium_Results.log"), skip = 20)
    LDSC_Consortium_Results <- as.data.frame(LDSC_Consortium_Results)
    LDSC_Consortium_Intercept <- gsub("Intercept: ","",LDSC_Consortium_Results[str_detect(LDSC_Consortium_Results[,1],"Intercept:"),1])
    LDSC_Consortium_h2 <- gsub("Total Observed scale h2: ","",LDSC_Consortium_Results[str_detect(LDSC_Consortium_Results[,1],"Total Observed scale h2:"),1])
  }else{
    LDSC_Consortium_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Consortium_Results.log"), skip = 20)
    LDSC_Consortium_Results <- as.data.frame(LDSC_Consortium_Results)
    LDSC_Consortium_Intercept <- gsub("Intercept: ","",LDSC_Consortium_Results[str_detect(LDSC_Consortium_Results[,1],"Intercept:"),1])
    LDSC_Consortium_h2 <- gsub("Total Observed scale h2: ","",LDSC_Consortium_Results[str_detect(LDSC_Consortium_Results[,1],"Total Observed scale h2:"),1]) 
  }
  
  if(trait %in% c("LDL","TC")){
    LDSC_Imputed_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_statin_adj_LDSC_Imputed_Results.log"), skip = 20)
    LDSC_Imputed_Results <- as.data.frame(LDSC_Imputed_Results)
    LDSC_Imputed_Intercept <- gsub("Intercept: ","",LDSC_Imputed_Results[str_detect(LDSC_Imputed_Results[,1],"Intercept:"),1])
    LDSC_Imputed_h2 <- gsub("Total Observed scale h2: ","",LDSC_Imputed_Results[str_detect(LDSC_Imputed_Results[,1],"Total Observed scale h2:"),1])
  }else{
    LDSC_Imputed_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Imputed_Results.log"), skip = 20)
    LDSC_Imputed_Results <- as.data.frame(LDSC_Imputed_Results)
    LDSC_Imputed_Intercept <- gsub("Intercept: ","",LDSC_Imputed_Results[str_detect(LDSC_Imputed_Results[,1],"Intercept:"),1])
    LDSC_Imputed_h2 <- gsub("Total Observed scale h2: ","",LDSC_Imputed_Results[str_detect(LDSC_Imputed_Results[,1],"Total Observed scale h2:"),1])
  }
  
  LDSC_NewPCs_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_newpcs_LDSC_Results.log"), skip = 20)
  LDSC_NewPCs_Results <- as.data.frame(LDSC_NewPCs_Results)
  LDSC_NewPCs_Intercept <- gsub("Intercept: ","",LDSC_NewPCs_Results[str_detect(LDSC_NewPCs_Results[,1],"Intercept:"),1])
  LDSC_NewPCs_h2 <- gsub("Total Observed scale h2: ","",LDSC_NewPCs_Results[str_detect(LDSC_NewPCs_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_RankNormal_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"adj_norm_ranknormal_LDSC_Results.log"), skip = 20)
  LDSC_RankNormal_Results <- as.data.frame(LDSC_RankNormal_Results)
  LDSC_RankNormal_Intercept <- gsub("Intercept: ","",LDSC_RankNormal_Results[str_detect(LDSC_RankNormal_Results[,1],"Intercept:"),1])
  LDSC_RankNormal_h2 <- gsub("Total Observed scale h2: ","",LDSC_RankNormal_Results[str_detect(LDSC_RankNormal_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Origianl_LD_Results <- read_csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Original_Results.log"), skip = 20)
  LDSC_Origianl_LD_Results <- as.data.frame(LDSC_Origianl_LD_Results)
  LDSC_Origianl_LD_Intercept <- gsub("Intercept: ","",LDSC_Origianl_LD_Results[str_detect(LDSC_Origianl_LD_Results[,1],"Intercept:"),1])
  LDSC_Origianl_LD_h2 <- gsub("Total Observed scale h2: ","",LDSC_Origianl_LD_Results[str_detect(LDSC_Origianl_LD_Results[,1],"Total Observed scale h2:"),1])
  
  LDSC_Results <- rbind(LDSC_Results,data.frame(trait = trait,method = c("plink 20 PCs","plink 10 PCs (original)","regenie 10 PCs","consortium","imputed 10 PCs (new data)","plink new 10 PCs","plink rank normal 10 PCs","plink 10 PCs original LD"),
                                                h2 = c(LDSC_20PCs_h2,LDSC_Original_h2,LDSC_Regenie_h2,LDSC_Consortium_h2,LDSC_Imputed_h2,LDSC_NewPCs_h2,LDSC_RankNormal_h2,LDSC_Origianl_LD_h2),
                                                Intercept = c(LDSC_20PCs_Intercept,LDSC_Original_Intercept,LDSC_Regenie_Intercept,LDSC_Consortium_Intercept,LDSC_Imputed_Intercept,LDSC_NewPCs_Intercept,LDSC_RankNormal_Intercept,LDSC_Origianl_LD_Intercept)))
  
}

write.csv(LDSC_Results,file = "LDSC_Results.csv",row.names = FALSE)