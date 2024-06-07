rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/Y_n_140488_h2_common_0.05_h2_rare_0.0125.RData")

set.seed(1335)

sampleids_all <- Y[[1]]$IDs

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids_all,]

# sum(ukb_pheno$IID == sampleids_all) == nrow(ukb_pheno); TRUE

i <- (1:length(sampleids_all))[ukb_pheno$ancestry == "EUR"]

train_number <- round(length(sampleids_all)*0.7) + 1
train <- sample(i, train_number)

i <- (1:length(sampleids_all))[!((1:length(sampleids_all)) %in% train)]
i_EUR <- i[ukb_pheno$ancestry[i] == "EUR"]
i_AFR <- i[ukb_pheno$ancestry[i] == "AFR"]
i_SAS <- i[ukb_pheno$ancestry[i] == "SAS"]
i_EAS <- i[ukb_pheno$ancestry[i] == "EAS"]
i_AMR <- i[ukb_pheno$ancestry[i] == "AMR"]

tune <- c(sample(i_EUR,round(length(i_EUR)/2)),
          sample(i_AFR,round(length(i_AFR)/2)),
          sample(i_SAS,round(length(i_SAS)/2)),
          sample(i_EAS,round(length(i_EAS)/2)),
          sample(i_AMR,round(length(i_AMR)/2)))

validation <- i[!(i %in% tune)]

train <- sampleids_all[train]
tune <- sampleids_all[tune]
validation <- sampleids_all[validation]

write.table(train,"/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/train.txt",row.names = FALSE,col.names = FALSE)
write.table(tune,"/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/tune.txt",row.names = FALSE,col.names = FALSE)
write.table(validation,"/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/validation.txt",row.names = FALSE,col.names = FALSE)

reference <- train[sample(1:length(train),3000,replace = FALSE)]

write.table(cbind(0,reference),"/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/reference_CT.txt",row.names = FALSE,col.names = FALSE)
write.table(reference,"/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/reference.txt",row.names = FALSE,col.names = FALSE)

Y_train <- list()
Y_tune <- list()
Y_validation <- list()

for(i in 1:length(Y)){
  Y_train[[i]] <- Y[[i]][Y[[i]]$IDs %in% train,]
  Y_tune[[i]] <- Y[[i]][Y[[i]]$IDs %in% tune,]
  Y_validation[[i]] <- Y[[i]][Y[[i]]$IDs %in% validation,]
}

save(Y_train,file = "/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")
save(Y_tune,file = "/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")
save(Y_validation,file = "/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")

system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/reference_CT.txt --make-bed --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/reference")

library(bigsnpr)
if(file.exists("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.rds")){
  file.remove("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.rds")
  file.remove("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.bk")
}

#### read in reference data, this should match as this is what the reference data was in CT
snp_readBed("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.bed",backingfile = "/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference")






rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/Y_n_140488_h2_common_0.05_h2_rare_0.0125.RData")

set.seed(1335)

sampleids_all <- Y[[1]]$IDs

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids_all,]

# sum(ukb_pheno$IID == sampleids_all) == nrow(ukb_pheno); TRUE

i <- (1:length(sampleids_all))[ukb_pheno$ancestry == "EUR"]

train_number <- round(length(sampleids_all)*0.35) + 1
train <- sample(i, train_number)

i <- (1:length(sampleids_all))[!((1:length(sampleids_all)) %in% train)]
i_EUR <- i[ukb_pheno$ancestry[i] == "EUR"]
i_AFR <- i[ukb_pheno$ancestry[i] == "AFR"]
i_SAS <- i[ukb_pheno$ancestry[i] == "SAS"]
i_EAS <- i[ukb_pheno$ancestry[i] == "EAS"]
i_AMR <- i[ukb_pheno$ancestry[i] == "AMR"]

tune <- c(sample(i_EUR,round(length(i_EUR)/2)),
          sample(i_AFR,round(length(i_AFR)/2)),
          sample(i_SAS,round(length(i_SAS)/2)),
          sample(i_EAS,round(length(i_EAS)/2)),
          sample(i_AMR,round(length(i_AMR)/2)))

validation <- i[!(i %in% tune)]

train <- sampleids_all[train]
tune <- sampleids_all[tune]
validation <- sampleids_all[validation]

write.table(train,"/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/train.txt",row.names = FALSE,col.names = FALSE)
write.table(tune,"/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/tune.txt",row.names = FALSE,col.names = FALSE)
write.table(validation,"/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/validation.txt",row.names = FALSE,col.names = FALSE)

reference <- train[sample(1:length(train),3000,replace = FALSE)]

write.table(cbind(0,reference),"/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/reference_CT.txt",row.names = FALSE,col.names = FALSE)
write.table(reference,"/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/reference.txt",row.names = FALSE,col.names = FALSE)

Y_train <- list()
Y_tune <- list()
Y_validation <- list()

for(i in 1:length(Y)){
  Y_train[[i]] <- Y[[i]][Y[[i]]$IDs %in% train,]
  Y_tune[[i]] <- Y[[i]][Y[[i]]$IDs %in% tune,]
  Y_validation[[i]] <- Y[[i]][Y[[i]]$IDs %in% validation,]
}

save(Y_train,file = "/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")
save(Y_tune,file = "/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")
save(Y_validation,file = "/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")

system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/reference_CT.txt --make-bed --out /data/williamsjacr/UKB_WES_Simulation/Simulation2/reference")

library(bigsnpr)
if(file.exists("/data/williamsjacr/UKB_WES_Simulation/Simulation2/reference.rds")){
  file.remove("/data/williamsjacr/UKB_WES_Simulation/Simulation2/reference.rds")
  file.remove("/data/williamsjacr/UKB_WES_Simulation/Simulation2/reference.bk")
}

#### read in reference data, this should match as this is what the reference data was in CT
snp_readBed("/data/williamsjacr/UKB_WES_Simulation/Simulation2/reference.bed",backingfile = "/data/williamsjacr/UKB_WES_Simulation/Simulation2/reference")

