rm(list = ls())

sampleids_all <- read.table("/data/williamsjacr/UKB_WES_lipids/Data/sampleids.txt", quote="\"", comment.char="")

unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")

load("/data/BB_Bioinformatics/ProjectData/UKB_WES_lipids/tmp.LDL.20211014.Rdata")
phenotype <- tmp.LDL
rm(tmp.LDL)

save(phenotype,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_All.RData")

sampleids_all <- sampleids_all[sampleids_all[,1] %in% phenotype$userId,,drop = FALSE]
sampleids_all <- sampleids_all[sampleids_all[,1] %in% unrels_nRandomSNPs_0$V1,,drop = FALSE]

set.seed(1335)

i <- 1:nrow(sampleids_all)

train_number <- round(max(i)*0.7) + 1
tune_val_number <- round(max(i)*0.15)

train <- sample(i, train_number)

i <- i[!(i %in% train)]

tune <- sample(i, tune_val_number)
validation <- i[!(i %in% tune)]

train <- sampleids_all[train,]
tune <- sampleids_all[tune,]
validation <- sampleids_all[validation,]

write.table(train,"/data/williamsjacr/UKB_WES_lipids/Data/train.txt",row.names = FALSE,col.names = FALSE)
write.table(tune,"/data/williamsjacr/UKB_WES_lipids/Data/tune.txt",row.names = FALSE,col.names = FALSE)
write.table(validation,"/data/williamsjacr/UKB_WES_lipids/Data/validation.txt",row.names = FALSE,col.names = FALSE)

reference <- train[sample(1:length(train),3000,replace = FALSE)]

write.table(reference,"/data/williamsjacr/UKB_WES_lipids/Data/reference.txt",row.names = FALSE,col.names = FALSE)

phenotype_train <- phenotype[phenotype$userId %in% train,]
phenotype_tune <- phenotype[phenotype$userId %in% tune,]
phenotype_validation <- phenotype[phenotype$userId %in% validation,]

phenotype_train$FID <- 0
phenotype_tune$FID <- 0
phenotype_validation$FID <- 0

phenotype_train <- phenotype_train[,c("userId","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
phenotype_tune <- phenotype_tune[,c("userId","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
phenotype_validation <- phenotype_validation[,c("userId","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]

colnames(phenotype_train) <- c("IID",colnames(phenotype_train)[-1])
colnames(phenotype_tune) <- c("IID",colnames(phenotype_tune)[-1])
colnames(phenotype_validation) <- c("IID",colnames(phenotype_validation)[-1])

phenotype_train[,c(4:5,7:16)] <- scale(phenotype_train[,c(4:5,7:16)]) 
phenotype_tune[,c(4:5,7:16)] <- scale(phenotype_tune[,c(4:5,7:16)]) 
phenotype_validation[,c(4:5,7:16)] <- scale(phenotype_validation[,c(4:5,7:16)]) 

save(phenotype_train,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.RData")
save(phenotype_tune,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.RData")
save(phenotype_validation,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.RData")

write.table(phenotype_train,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt",sep = '\t',row.names = FALSE,quote = FALSE)
write.table(phenotype_tune,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt",sep = '\t',row.names = FALSE,quote = FALSE)
write.table(phenotype_validation,file = "/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt",sep = '\t',row.names = FALSE,quote = FALSE)

