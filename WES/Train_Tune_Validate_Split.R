rm(list = ls())

library(bigsnpr)

sampleids_all <- read.table("/data/williamsjacr/UKB_WES_Full_Processed_Data/sampleids.txt", quote="\"", comment.char="")

unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

ukb_pheno <- as.data.frame(ukb_pheno)

ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids_all[,1],]

sampleids_all <- sampleids_all[sampleids_all[,1] %in% ukb_pheno$IID,,drop = FALSE]
sampleids_all <- sampleids_all[sampleids_all[,1] %in% unrels_nRandomSNPs_0$V1,,drop = FALSE]

set.seed(1335)

i <- (1:nrow(sampleids_all))[ukb_pheno$ancestry == "EUR"]

train_number <- round(length(i)*0.7) + 1
tune_val_number <- round(length(i)*0.15)

train <- sample(i, train_number)

i <- i[!(i %in% train)]

tune <- sample(i, tune_val_number)
validation <- c(i[!(i %in% tune)],(1:length(sampleids_all))[ukb_pheno$ancestry != "EUR"])

train <- sampleids_all[train,]
tune <- sampleids_all[tune,]
validation <- sampleids_all[validation,]

write.table(train,"/data/williamsjacr/UKB_WES_Phenotypes/train.txt",row.names = FALSE,col.names = FALSE)
write.table(tune,"/data/williamsjacr/UKB_WES_Phenotypes/tune.txt",row.names = FALSE,col.names = FALSE)
write.table(validation,"/data/williamsjacr/UKB_WES_Phenotypes/validation.txt",row.names = FALSE,col.names = FALSE)

reference <- train[sample(1:length(train),3000,replace = FALSE)]

write.table(reference,"/data/williamsjacr/UKB_WES_Phenotypes/reference.txt",row.names = FALSE,col.names = FALSE)

for(i in 1:22){
  system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr",i,"/ukbb_wes_200k_chr",i,"_common --keep /data/williamsjacr/UKB_WES_Phenotypes/reference.txt --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/ukbb_wes_200k_chr",i,"_common_reference"))
}

system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/reference.txt --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_reference"))
# system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_train"))
# system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_tune"))
# system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_validation"))


phenotype_train <- ukb_pheno[ukb_pheno$IID %in% train,]
phenotype_tune <- ukb_pheno[ukb_pheno$IID %in% tune,]
phenotype_validation <- ukb_pheno[ukb_pheno$IID %in% validation,]

phenotype_train$FID <- 0
phenotype_tune$FID <- 0
phenotype_validation$FID <- 0

phenotype_train <- phenotype_train[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
phenotype_train$age2 <- phenotype_train$age^2
phenotype_tune <- phenotype_tune[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
phenotype_tune$age2 <- phenotype_tune$age^2
phenotype_validation <- phenotype_validation[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
phenotype_validation$age2 <- phenotype_validation$age^2

colnames(phenotype_train) <- c("IID",colnames(phenotype_train)[-1])
colnames(phenotype_tune) <- c("IID",colnames(phenotype_tune)[-1])
colnames(phenotype_validation) <- c("IID",colnames(phenotype_validation)[-1])

phenotype_train[c(14,16:26)] <- lapply(phenotype_train[c(14,16:26)], function(x){c(scale(x))})
phenotype_tune[c(14,16:26)] <- lapply(phenotype_tune[c(14,16:26)], function(x){c(scale(x))})
phenotype_validation[c(14,16:26)] <- lapply(phenotype_validation[c(14,16:26)], function(x){c(scale(x))})

save(phenotype_train,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Train.RData")
save(phenotype_tune,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.RData")
save(phenotype_validation,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.RData")

write.table(phenotype_train,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt",sep = '\t',row.names = FALSE,quote = FALSE)
write.table(phenotype_tune,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt",sep = '\t',row.names = FALSE,quote = FALSE)
write.table(phenotype_validation,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt",sep = '\t',row.names = FALSE,quote = FALSE)

snp_readBed("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_reference.bed",backingfile = "/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_reference")
