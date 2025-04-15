rm(list=ls())

library(dplyr)
library(stringr)

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")

binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")

trait <- as.numeric(commandArgs(TRUE)[1])

trait <- c(continuous_traits,binary_traits)[trait]

if(trait %in% continuous_traits){
  dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_",trait,".regenie"), sep="")
  colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
  dat$P <- 10^(-1*dat$LOG10P)
  
  dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
  colnames(dat) <- c("CHR","SNP","REF","BP","ALT","BETA","P")
}else{
  if(trait %in% c("Breast","Prostate")){
    fill <- "bp"
    confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
  }else{
    fill <- "act"
    confounders <-  paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
  }
  
  dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_",fill,"_",trait,".regenie"), sep="")
  colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
  dat$P <- 10^(-1*dat$LOG10P)
  
  dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
  colnames(dat) <- c("CHR","SNP","REF","BP","ALT","BETA","P")
}

beta_final <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_CT_SNP_List.csv"))

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
dat$CT_p_value_1 <- 0
dat$CT_p_value_1[dat$P <= pthres[1]] <- dat$BETA[dat$P <= pthres[1]]
dat$CT_p_value_1[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_2 <- 0
dat$CT_p_value_2[dat$P <= pthres[2]] <- dat$BETA[dat$P <= pthres[2]]
dat$CT_p_value_2[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_3 <- 0
dat$CT_p_value_3[dat$P <= pthres[3]] <- dat$BETA[dat$P <= pthres[3]]
dat$CT_p_value_3[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_4 <- 0
dat$CT_p_value_4[dat$P <= pthres[4]] <- dat$BETA[dat$P <= pthres[4]]
dat$CT_p_value_4[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_5 <- 0
dat$CT_p_value_5[dat$P <= pthres[5]] <- dat$BETA[dat$P <= pthres[5]]
dat$CT_p_value_5[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_6 <- 0
dat$CT_p_value_6[dat$P <= pthres[6]] <- dat$BETA[dat$P <= pthres[6]]
dat$CT_p_value_6[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_7 <- 0
dat$CT_p_value_7[dat$P <= pthres[7]] <- dat$BETA[dat$P <= pthres[7]]
dat$CT_p_value_7[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_8 <- 0
dat$CT_p_value_8[dat$P <= pthres[8]] <- dat$BETA[dat$P <= pthres[8]]
dat$CT_p_value_8[!(dat$SNP %in% beta_final$SNP)] <- 0
dat$CT_p_value_9 <- 0
dat$CT_p_value_9[dat$P <= pthres[9]] <- dat$BETA[dat$P <= pthres[9]]
dat$CT_p_value_9[!(dat$SNP %in% beta_final$SNP)] <- 0

Final_Coefficients_CT <- dat[,c("CHR","SNP","REF","BP","ALT","P",paste0("CT_p_value_",1:9))]
prs.file <- Final_Coefficients_CT[,c("SNP","ALT",paste0("CT_p_value_",1:9))]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/test_validation_",trait))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
prs_all_validation <- prs_all_validation[,-1]
colnames(test_validation) <- colnames(prs_all_validation)
all.equal(test_validation,prs_all_validation)
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/test_validation_",trait,".sscore"))
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Final_Score"))

Final_Coefficients_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt"), sep="")
Final_Coefficients_LDPred <- left_join(Final_Coefficients_CT[,c("SNP","ALT")],Final_Coefficients_LDPred)
Final_Coefficients_LDPred[is.na(Final_Coefficients_LDPred)] <- 0

## Works
prs.file <- Final_Coefficients_LDPred[,c("SNP","ALT",paste0("BETA.e",1:255))]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/test_validation_",trait))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation.sscore"), sep="")
prs_all_validation <- prs_all_validation[,c(2,5:ncol(prs_all_validation))]
colnames(test_validation) <- colnames(prs_all_validation)
all.equal(test_validation,prs_all_validation)
cor(test_validation$SCORE7_SUM,prs_all_validation$SCORE7_SUM)
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/test_validation_",trait,".sscore"))
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Final_Score"))

Final_Coefficients_LASSOSUM <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt"), sep="")
Final_Coefficients_LASSOSUM <- left_join(Final_Coefficients_CT[,c("SNP","ALT")],Final_Coefficients_LASSOSUM)
Final_Coefficients_LASSOSUM[is.na(Final_Coefficients_LASSOSUM)] <- 0
prs.file <- Final_Coefficients_LASSOSUM[,c("SNP","ALT",paste0("BETA.e",1:300))]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/test_validation_",trait))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation.sscore"), sep="")
prs_all_validation <- prs_all_validation[,c(2,5:ncol(prs_all_validation))]
colnames(test_validation) <- colnames(prs_all_validation)
all.equal(test_validation,prs_all_validation)
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/test_validation_",trait,".sscore"))
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Final_Score"))

#### Linear Model Way
CT_prs_all_tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
LDpred2_prs_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune.sscore"), header=FALSE, comment.char="#")
LASSOsum2_prs_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune.sscore"), header=FALSE, comment.char="#")

prs_tune_all <- cbind(LDpred2_prs_tune[,1:2],CT_prs_all_tune[,-c(1,2)],LDpred2_prs_tune[,-c(1,2,3,4)],LASSOsum2_prs_tune[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("FID","IID",paste0("PRS",1:(ncol(prs_tune_all) - 2),"_",trait))
rm(CT_prs_all_tune);rm(LDpred2_prs_tune);rm(LASSOsum2_prs_tune)

## Merge covariates and y for tuning with the prs_mat
Best_Tune_All <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Tune.csv"))

# PRS_Tune <- as.matrix(pheno_tune[,names(Final_Coefficients_Ensemble$Coefficients)[-1]]) %*% matrix(Final_Coefficients_Ensemble$Coefficients[-1],ncol = 1)
# cor(PRS_Tune,Best_Tune_All$PRS)

Best_Tune_All <- inner_join(Best_Tune_All,prs_tune_all)

Beta_Star <- matrix(unname(coef(lm(PRS ~.,data = Best_Tune_All[,c("PRS",colnames(Best_Tune_All)[str_detect(colnames(Best_Tune_All),paste0("_",trait))])]))),ncol = 1)
Beta_Star <- Beta_Star[-1,,drop = FALSE]
Beta_Star[is.na(Beta_Star)] <- 0

All_Scores <- cbind(Final_Coefficients_CT[,7:ncol(Final_Coefficients_CT)],Final_Coefficients_LDPred[,4:ncol(Final_Coefficients_LDPred)],Final_Coefficients_LASSOSUM[,4:ncol(Final_Coefficients_LASSOSUM)])

Score_Final <- as.matrix(All_Scores) %*% Beta_Star

prs.file <- data.frame(SNP = Final_Coefficients_CT$SNP,ALT = Final_Coefficients_CT$ALT,BETA = Score_Final)
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/test_validation_",trait))

Best_Validation_All <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
test_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/test_validation_",trait,".sscore"), header=FALSE, comment.char="#")

all.equal(Best_Validation_All$IID,test_validation[,2])
all.equal(Best_Validation_All$PRS,test_validation[,5])
cor(Best_Validation_All$PRS,test_validation[,5])

system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/test_validation_",trait,".sscore"))
system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/test_validation_",trait,".log"))

