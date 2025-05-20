rm(list=ls())

library(dplyr)
library(stringr)

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "Asthma"
}else if(trait == 2){
  trait <- "CAD"
}else if(trait == 3){
  trait <- "T2D"
}else if(trait == 4){
  trait <- "Breast"
}else{
  trait <- "Prostate"
}

Final_Coefficients_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_Final_Coefficients.csv"))
prs.file <- Final_Coefficients_CT[,c("SNP","A1",paste0("CT_p_value_",1:9))]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
prs_all_validation <- prs_all_validation[,-1]
colnames(test_validation) <- colnames(prs_all_validation)
all.equal(test_validation,prs_all_validation)

Final_Coefficients_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_Final_Coefficients.csv"))
prs.file <- Final_Coefficients_LDPred[,c("SNP","A1",paste0("LDPred2_SCORE",1:255,"_SUM"))]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_validation.sscore"), sep="")
prs_all_validation <- prs_all_validation[,c(2,5:ncol(prs_all_validation))]
colnames(test_validation) <- colnames(prs_all_validation)
all.equal(test_validation,prs_all_validation)

Final_Coefficients_LASSOSUM <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_Final_Coefficients.csv"))
prs.file <- Final_Coefficients_LASSOSUM[,c("SNP","A1",paste0("LASSOSum2_SCORE",1:300,"_SUM"))]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_validation.sscore"), sep="")
prs_all_validation <- prs_all_validation[,c(2,5:ncol(prs_all_validation))]
colnames(test_validation) <- colnames(prs_all_validation)
all.equal(test_validation,prs_all_validation)

#### Linear Model Way
prs_tune_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
prs_tune_LDPred2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_tune.sscore"))
prs_tune_LASSOSum2 <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_tune.sscore"))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
Best_Tune_All <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Tune_All.txt"))
Best_Tune_All <- left_join(Best_Tune_All,prs_tune_all,by = "IID")
Best_Tune_All <- subset(Best_Tune_All,select = -c(IID,X.FID))

Beta_Star <- matrix(unname(coef(lm(prs ~.,data = Best_Tune_All))),ncol = 1)
Beta_Star <- Beta_Star[-1,,drop = FALSE]
Beta_Star[is.na(Beta_Star),1] <- 0

All <- cbind(Final_Coefficients_CT,Final_Coefficients_LDPred[,str_detect(colnames(Final_Coefficients_LDPred),"LDPred2")],Final_Coefficients_LASSOSUM[,str_detect(colnames(Final_Coefficients_LASSOSUM),"LASSOSum2")])
score_full <- All[,colnames(All) %in% names(coef(lm(prs ~.,data = Best_Tune_All)))[-1]]
score_full <- score_full[,match(colnames(score_full),names(coef(lm(prs ~.,data = Best_Tune_All)))[-1])]
score_full <- as.matrix(score_full)

prs.file <- data.frame(SNP = All$SNP,ALT = All$A1,BETA = score_full %*% Beta_Star)
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait))

Best_Validation_All <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")

all.equal(Best_Validation_All$IID,test_validation[,2])
all.equal(Best_Validation_All$prs,test_validation[,5])
cor(Best_Validation_All$prs,test_validation[,5])

system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait,".sscore"))
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Binary/test_validation_",trait,".log"))

