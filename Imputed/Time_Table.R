rm(list = ls())

dat <- read.table(textConnection("start.date start.time end.date end.time
2025-03-26   09:38:56 2025-03-26   17:53:46"), header=TRUE)
regenie_time_continuous <-  as.numeric(difftime(strptime(paste(dat[,3],dat[,4]),"%Y-%m-%d %H:%M:%S"),strptime(paste(dat[,1],dat[,2]),"%Y-%m-%d %H:%M:%S"),units = "hours"))

dat <- read.table(textConnection("start.date start.time end.date end.time
2025-03-27   21:29:52 2025-03-29   01:11:37 
2025-03-29   01:11:38 2025-03-29   03:28:15"), header=TRUE)
regenie_time_binary_withsex <-  sum(as.numeric(difftime(strptime(paste(dat[,3],dat[,4]),"%Y-%m-%d %H:%M:%S"),strptime(paste(dat[,1],dat[,2]),"%Y-%m-%d %H:%M:%S"),units = "hours")))

dat <- read.table(textConnection("start.date start.time end.date end.time
2025-03-27   21:29:52 2025-03-28   23:18:20
2025-03-29   01:11:38 2025-03-29   03:23:58"), header=TRUE)
regenie_time_binary_withoutsex <-  sum(as.numeric(difftime(strptime(paste(dat[,3],dat[,4]),"%Y-%m-%d %H:%M:%S"),strptime(paste(dat[,1],dat[,2]),"%Y-%m-%d %H:%M:%S"),units = "hours")))

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

Time_Data <- data.frame(Trait = c("BMI","Height","LDL","HDL","logTG","TC","Asthma","Breast","CAD","Prostate","T2D"),GWAS = c(rep(regenie_time_continuous,6),c(regenie_time_binary_withsex,regenie_time_binary_withoutsex,regenie_time_binary_withsex,regenie_time_binary_withoutsex,regenie_time_binary_withsex)), CT = CT_Times, LDpred2_Lassosum2 = LDpred2_Lassosum2_Times,RICE_CV = RICE_CV_Times,STAARpipeline = STAARpipeline_Times, RICE_RV = RICE_RV_Times,RICE = RICE_Times)
Time_Data$Total <- rowSums(Time_Data[,-1])

for(i in colnames(Time_Data)[-c(1,ncol(Time_Data))]){
  Time_Data[,i] <- paste0(round(Time_Data[,i],digits = 2),"(",round(100*Time_Data[,i]/Time_Data$Total,digits = 2),"%)")
}
Time_Data$Total <- as.character(round(Time_Data$Total,digits = 2))

write.csv(Time_Data,file = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Timings.csv",row.names = FALSE)
