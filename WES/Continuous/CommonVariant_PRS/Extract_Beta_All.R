rm(list=ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)
library(stringr)

BMI_Final_Coefficients_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/","BMI","_Final_Coefficients.csv"))
prs.file <- BMI_Final_Coefficients_CT[,c("SNP","A1",paste0("CT_p_value_",1:9))]
write.table(prs.file,file = "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score",col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out test_validation"))
test_validation <- read.delim("test_validation.sscore", header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
BMI_prs_all_validation <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/BMI_prs_all_validation.txt", sep="")
BMI_prs_all_validation <- BMI_prs_all_validation[,-1]
all.equal(test_validation,BMI_prs_all_validation)

BMI_Final_Coefficients_LDPred <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/BMI_Final_Coefficients.csv")
prs.file <- BMI_Final_Coefficients_LDPred[,c("SNP","A1",paste0("LDPred2_SCORE",1:255,"_SUM"))]
write.table(prs.file,file = "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score",col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out test_validation"))
test_validation <- read.delim("test_validation.sscore", header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
BMI_prs_all_validation <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/BMI_prs_validation.sscore", sep="")
BMI_prs_all_validation <- BMI_prs_all_validation[,c(2,5:ncol(BMI_prs_all_validation))]
colnames(test_validation) <- colnames(BMI_prs_all_validation)
all.equal(test_validation,BMI_prs_all_validation)

BMI_Final_Coefficients_LASSOSUM <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/BMI_Final_Coefficients.csv")
prs.file <- BMI_Final_Coefficients_LASSOSUM[,c("SNP","A1",paste0("LASSOSum2_SCORE",1:300,"_SUM"))]
write.table(prs.file,file = "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score",col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out test_validation"))
test_validation <- read.delim("test_validation.sscore", header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
BMI_prs_all_validation <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/BMI_prs_validation.sscore", sep="")
BMI_prs_all_validation <- BMI_prs_all_validation[,c(2,5:ncol(BMI_prs_all_validation))]
colnames(test_validation) <- colnames(BMI_prs_all_validation)
all.equal(test_validation,BMI_prs_all_validation)


BMI_All <- cbind(BMI_Final_Coefficients_CT,BMI_Final_Coefficients_LDPred[,str_detect(colnames(BMI_Final_Coefficients_LDPred),"LDPred2")],BMI_Final_Coefficients_LASSOSUM[,str_detect(colnames(BMI_Final_Coefficients_LASSOSUM),"LASSOSum2")])

BMI_Final_Coefficients_SL <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/BMI_final_coef.csv")

CV_SL_Beta <- BMI_Final_Coefficients_CT[,c("CHR","SNP","REF","BP","A1")]
CV_SL_Beta$BETA <- 0

for(i in 1:nrow(BMI_Final_Coefficients_SL)){
  name_i <- BMI_Final_Coefficients_SL$Coef[i]
  if(name_i %in% colnames(BMI_All)){
    CV_SL_Beta$BETA <- CV_SL_Beta$BETA + BMI_All[,name_i]*BMI_Final_Coefficients_SL$Beta[i] 
  }
}

prs.file <- CV_SL_Beta[,c("SNP","A1","BETA")]
write.table(prs.file,file = "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score",col.names = T,row.names = F,quote=F)

system("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Continuous/BMI_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out test_validation")

Best_Validation_All <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/BMI_Best_Validation_All.txt")
test_validation <- read.delim("test_validation.sscore", header=FALSE, comment.char="#")

all.equal(Best_Validation_All$IID,test_validation[,2])
all.equal(Best_Validation_All$prs,test_validation[,5])
cor(Best_Validation_All$prs,test_validation[,5])


BMI_Final_Coefficients_RV <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/BMI_final_coef.csv")






























## CT

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
names_CT <- paste0("p_value_",1:length(pthres))
CT_allcoef <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/all_chr_assoc.txt", sep="")
CT_allcoef <- CT_allcoef[,c("CHR","SNP","BP","REF","A1","BETA","P")]
colnames(CT_allcoef)[7] <- "P_Original"
CT_clump_pvals <- as.data.frame(fread("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/CT/LD_clump.clumped"))
CT_clump_pvals <- CT_clump_pvals[,c("CHR","SNP","BP","P")]

CT_something <- left_join(CT_allcoef,CT_clump_pvals)

CT_something$P_Original[is.na(CT_something$P)] <- NA

for(i in 1:length(pthres)){
  CT_something[paste0("CT_",names_CT[i])] <- ifelse(CT_something$P_Original <= pthres[i],CT_something$BETA,0)
  CT_something[paste0("CT_",names_CT[i])][is.na(CT_something[paste0("CT_",names_CT[i])])] <- 0
}

CT_something <- CT_something[,c("CHR","SNP","BP","REF","A1",paste0("CT_",names_CT))]

rm(list=setdiff(ls(), "CT_something"))

## LDPred

LDPred_something <- NULL

for(i in 1:22){
  LDPred_something <- rbind(LDPred_something,read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2-chr",i,".txt"), sep=""))
}

CT_allcoef <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/all_chr_assoc.txt", sep="")
CT_allcoef <- CT_allcoef[,c("CHR","SNP","BP","REF","A1")]

colnames(LDPred_something)[2:3] <- c("LDPred_A1","LDPred_REF")

LDPred_something <- left_join(CT_allcoef,LDPred_something)

LDPred_something$LDPred_A1[is.na(LDPred_something$LDPred_A1)] <- LDPred_something$A1[is.na(LDPred_something$LDPred_A1)]
LDPred_something$LDPred_REF[is.na(LDPred_something$LDPred_REF)] <- LDPred_something$REF[is.na(LDPred_something$LDPred_REF)]

# sum(LDPred_something$REF != LDPred_something$LDPred_REF)
# list1 <- which(LDPred_something$REF != LDPred_something$LDPred_REF)
# sum(LDPred_something$A1 != LDPred_something$LDPred_A1)
# list2 <- which(LDPred_something$A1 != LDPred_something$LDPred_A1)
# sum(LDPred_something$A1 == LDPred_something$LDPred_REF)
# list3 <- which(LDPred_something$A1 == LDPred_something$LDPred_REF)
# sum(LDPred_something$REF == LDPred_something$LDPred_A1)
# list4 <- which(LDPred_something$REF == LDPred_something$LDPred_A1)
# 
# all.equal(list1,list2)
# all.equal(list1,list3)
# all.equal(list1,list4)

LDPred_something[,8:58][is.na(LDPred_something[,8:58])] <- 0

LDPred_something[which(LDPred_something$REF != LDPred_something$LDPred_REF),8:58] <- -1*LDPred_something[which(LDPred_something$REF != LDPred_something$LDPred_REF),8:58]

colnames(LDPred_something)[8:58] <- paste0("LDPred2_SCORE",1:51,"_SUM")

rm(list=setdiff(ls(), c("CT_something","LDPred_something")))

## LASSOSum

LASSOSum_something <- NULL

for(i in 1:22){
  LASSOSum_something <- rbind(LASSOSum_something,read.csv(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2-chr",i,".txt"), sep=""))
}

CT_allcoef <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/all_chr_assoc.txt", sep="")
CT_allcoef <- CT_allcoef[,c("CHR","SNP","BP","REF","A1")]

colnames(LASSOSum_something)[2:3] <- c("LASSOSum_A1","LASSOSum_REF")

LASSOSum_something <- left_join(CT_allcoef,LASSOSum_something)

LASSOSum_something$LASSOSum_A1[is.na(LASSOSum_something$LASSOSum_A1)] <- LASSOSum_something$A1[is.na(LASSOSum_something$LASSOSum_A1)]
LASSOSum_something$LASSOSum_REF[is.na(LASSOSum_something$LASSOSum_REF)] <- LASSOSum_something$REF[is.na(LASSOSum_something$LASSOSum_REF)]

LASSOSum_something[,8:307][is.na(LASSOSum_something[,8:307])] <- 0

LASSOSum_something[which(LASSOSum_something$REF != LASSOSum_something$LASSOSum_REF),8:307] <- -1*LASSOSum_something[which(LASSOSum_something$REF != LASSOSum_something$LASSOSum_REF),8:307]

colnames(LASSOSum_something)[8:307] <- paste0("LASSOSum2_SCORE",1:300,"_SUM")

rm(list=setdiff(ls(), c("CT_something","LDPred_something","LASSOSum_something")))

betas_all <- cbind(CT_something,LDPred_something[,c(8:58)],LASSOSum_something[,8:307])

rm(list=setdiff(ls(), c("betas_all")))

#############################################################################

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/final_coefs_All.RData")

betas <- as.matrix(betas_all[,colnames(betas_all)[colnames(betas_all) %in% names(final_coefs)]])
final_coefs <- matrix(final_coefs,ncol = 1)

beta_final <- betas %*% final_coefs

beta_final <- data.frame(CHR = betas_all$CHR,SNP = betas_all$SNP,BP = betas_all$BP,REF = betas_all$REF, A1 = betas_all$A1,BETA = beta_final)

load("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Common_plus_RareVariants/Coefficients_STAARO.RData")

beta_final$BETA <- beta_final$BETA*effects[1]

write.table(beta_final,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_All_Beta.txt",sep = "\t",row.names = FALSE)

rm(list=setdiff(ls(), c("beta_final")))

#############################################################################
# Validate

prs.file <- beta_final[,c("SNP","A1","BETA")]
write.table(prs.file,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_All_Beta_Subset",col.names = T,row.names = F,quote=F)

system("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_All_Beta_Subset cols=+scoresums,-scoreavgs header no-mean-imputation --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_validation --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/test_validation")

Best_Validation_All <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt")
test_validation <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/test_validation.sscore", header=FALSE, comment.char="#")

all.equal(Best_Validation_All$IID,test_validation[,2])
all.equal(Best_Validation_All$prs,test_validation[,5])
cor(Best_Validation_All$prs,test_validation[,5])


