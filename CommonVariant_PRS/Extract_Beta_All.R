rm(list=ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)

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


