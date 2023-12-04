rm(list = ls())
load("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.RData")

phenotype_train <- phenotype_train[,c(2,1,3:ncol(phenotype_train))]
write.table(phenotype_train,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt",sep = '\t',row.names = FALSE,quote = FALSE)