rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")

for(i in 1:length(Y_train)){
  Y_train[[i]]$FID <- Y_train[[i]]$IDs
  Y_train[[i]] <- Y_train[[i]][,c(1,3,2)]
  colnames(Y_train[[i]]) <- c("IID","FID","Y")
  write.table(Y_train[[i]],file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train",i,".txt"),sep = '\t',row.names = FALSE,quote = FALSE)
}

chr22_filtered_common <- read.delim("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam", header=FALSE)
chr22_filtered_common$V1 <- 0
write.table(chr22_filtered_common,col.names = FALSE,row.names = FALSE,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam")

for(i in 1:length(Y_train)){
  system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --pheno /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train",i,".txt --pheno-name Y --linear allow-no-covars --vif 999 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/GWAS_Summary_Statistics/Y_Train",i))
}
