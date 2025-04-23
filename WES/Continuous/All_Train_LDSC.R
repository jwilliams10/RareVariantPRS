rm(list = ls())

library(bigsnpr)

sampleids_all <- read.table("/data/williamsjacr/UKB_WES_Full_Processed_Data/sampleids.txt", quote="\"", comment.char="")

unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

ukb_pheno <- as.data.frame(ukb_pheno)

ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids_all[,1],]

sampleids_all <- sampleids_all[sampleids_all[,1] %in% ukb_pheno$IID,,drop = FALSE]
sampleids_all <- sampleids_all[sampleids_all[,1] %in% unrels_nRandomSNPs_0$V1,,drop = FALSE]

sampleids_all <- sampleids_all[order(sampleids_all[,1]),,drop = FALSE]
ukb_pheno <- ukb_pheno[order(ukb_pheno$IID),]

set.seed(1330)

i <- (1:nrow(sampleids_all))[ukb_pheno$ancestry == "EUR"]

train_number <- round(nrow(sampleids_all)*0.7) + 1
train <- sample(i, train_number)

train <- sampleids_all[train,]

phenotype_train <- ukb_pheno[ukb_pheno$IID %in% train,]

phenotype_train$FID <- 0

phenotype_train <- phenotype_train[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex",paste0("pc",1:40))]
phenotype_train$age2 <- phenotype_train$age^2

colnames(phenotype_train) <- c("IID",colnames(phenotype_train)[-1])

phenotype_train[c(14,16:56)] <- lapply(phenotype_train[c(14,16:56)], function(x){c(scale(x))})

BMIadj_resid <- resid(lm(BMI~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = phenotype_train))
phenotype_train$BMIadj_norm[!is.na(phenotype_train$BMI)] <- qnorm((rank(BMIadj_resid,na.last="keep")-0.5)/length(BMIadj_resid))

Heightadj_resid <- resid(lm(Height~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = phenotype_train))
phenotype_train$Heightadj_norm[!is.na(phenotype_train$Height)] <- qnorm((rank(Heightadj_resid,na.last="keep")-0.5)/length(Heightadj_resid))

LDLadj_resid <- resid(lm(LDL~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = phenotype_train))
phenotype_train$LDLadj_norm[!is.na(phenotype_train$LDL)] <- qnorm((rank(LDLadj_resid,na.last="keep")-0.5)/length(LDLadj_resid))

HDLadj_resid <- resid(lm(HDL~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = phenotype_train))
phenotype_train$HDLadj_norm[!is.na(phenotype_train$HDL)] <- qnorm((rank(HDLadj_resid,na.last="keep")-0.5)/length(HDLadj_resid))

logTGadj_resid <- resid(lm(logTG~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = phenotype_train))
phenotype_train$logTGadj_norm[!is.na(phenotype_train$logTG)] <- qnorm((rank(logTGadj_resid,na.last="keep")-0.5)/length(logTGadj_resid))

TCadj_resid <- resid(lm(TC~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = phenotype_train))
phenotype_train$TCadj_norm[!is.na(phenotype_train$TC)] <- qnorm((rank(TCadj_resid,na.last="keep")-0.5)/length(TCadj_resid))

write.table(phenotype_train,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt",sep = '\t',row.names = FALSE,quote = FALSE)