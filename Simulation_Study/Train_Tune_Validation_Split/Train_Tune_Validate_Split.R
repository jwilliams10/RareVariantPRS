rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/Y_n_130368_h2_common_0.05_h2_rare_0.0125.RData")

set.seed(1335)

sampleids_all <- Y[[1]]$IDs

i <- 1:length(sampleids_all)

train_number <- round(max(i)*0.7) + 1
tune_val_number <- round(max(i)*0.15)

train <- sample(i, train_number)

i <- i[!(i %in% train)]

tune <- sample(i, tune_val_number)
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

