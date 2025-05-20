rm(list = ls())

dat <- read.table(textConnection("start.date start.time end.date end.time
2025-04-28 07:31:02 2025-04-28 16:03:51
2025-04-26 11:01:53 2025-04-26 18:33:09
2025-04-26 11:01:53 2025-04-26 23:40:08
2025-04-26 11:01:53 2025-04-26 20:53:02
2025-04-26 11:01:53 2025-04-26 19:56:16
2025-04-26 11:01:53 2025-04-26 18:04:33
2025-04-26 11:09:39 2025-04-26 23:37:28
2025-04-26 11:09:50 2025-04-26 21:56:40
2025-04-26 11:09:39 2025-04-26 22:37:45
2025-04-26 11:09:39 2025-04-26 22:50:14
2025-04-26 11:09:40 2025-04-26 23:53:01"), header=TRUE)
regenie_time <-  as.numeric(difftime(strptime(paste(dat[,3],dat[,4]),"%Y-%m-%d %H:%M:%S"),strptime(paste(dat[,1],dat[,2]),"%Y-%m-%d %H:%M:%S"),units = "hours"))

count <- 1
CT_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  CT_Times[count] <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",i,"_Time.RData")))/3600
  count <- count + 1
}

count <- 1
LDpred2_Lassosum2_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  LDpred2_Lassosum2_Times[count] <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/",i,"_LDpred2_Lassosum2_Time.RData")))/3600
  count <- count + 1
}

count <- 1
RICE_CV_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  RICE_CV_Times[count] <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",i,"_Time.RData")))/3600
  count <- count + 1
}

count <- 1
RICE_CV_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  RICE_CV_Times[count] <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",i,"_Time.RData")))/3600
  count <- count + 1
}

count <- 1
STAARpipeline_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  Null_Model_Time <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/NullModel_Time_",i,".RData")))/3600
  Analysis_Time <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",i,"_Analysis_Time.RData")))/3600
  Summary_Time <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",i,"_Summary_Time.RData")))/3600
  STAARpipeline_Times[count] <- Null_Model_Time + Analysis_Time + Summary_Time
  count <- count + 1
}

count <- 1
RICE_RV_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  RICE_RV_Times[count] <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",i,"_Time.RData")))/3600
  count <- count + 1
}

count <- 1
RICE_Times <- vector()
for(i in c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D")){
  RICE_Times[count] <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/",i,"_Time.RData")))/3600
  count <- count + 1
}

Time_Data <- data.frame(Trait = c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D"),GWAS = regenie_time, CT = CT_Times, LDpred2_Lassosum2 = LDpred2_Lassosum2_Times,RICE_CV = RICE_CV_Times,STAARpipeline = STAARpipeline_Times, RICE_RV = RICE_RV_Times,RICE = RICE_Times)
Time_Data$Total <- rowSums(Time_Data[,-1])

for(i in colnames(Time_Data)[-c(1,ncol(Time_Data))]){
  Time_Data[,i] <- paste0(round(Time_Data[,i],digits = 2),"(",round(100*Time_Data[,i]/Time_Data$Total,digits = 2),"%)")
}
Time_Data$Total <- as.character(round(Time_Data$Total,digits = 2))

write.csv(Time_Data,file = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Timings.csv",row.names = FALSE)
