rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train.RData")

i <- as.numeric(commandArgs(TRUE)[1])

Y_train[[i]]$FID <- Y_train[[i]]$IDs
Y_train[[i]] <- Y_train[[i]][,c(1,3,2)]
colnames(Y_train[[i]]) <- c("IID","FID","Y")
write.table(Y_train[[i]],file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train",i,".txt"),sep = '\t',row.names = FALSE,quote = FALSE)

system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --pheno /data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train",i,".txt --pheno-name Y --linear allow-no-covars --vif 999 --out /data/williamsjacr/UKB_WES_Simulation/Simulation4/GWAS_Summary_Statistics/Y_Train",i))

