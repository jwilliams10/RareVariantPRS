rm(list=ls())

library(dplyr)
library(stringr)

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/Extract_Betas_All_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/Extract_Betas_All_Binary.sh -icmd="bash Extract_Betas_All_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/ --priority low --instance-type mem3_ssd1_v2_x32
# done

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

Final_Coefficients_CT <- read.csv(paste0(trait,"_Final_Coefficients.csv"))
system(paste0("rm ",paste0(trait,"_Final_Coefficients.csv")))
prs.file <- Final_Coefficients_CT[,c("SNP","A1",paste0("CT_p_value_",1:9))]
write.table(prs.file,file = paste0(trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("plink2 --score ",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile all_chr --keep validation.txt --threads 1 --out test_validation_",trait))
test_validation <- read.delim(paste0("test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0(trait,"_prs_all_validation.txt"), sep="")
system(paste0("rm ",paste0(trait,"_prs_all_validation.txt")))
prs_all_validation <- prs_all_validation[,-1]
colnames(test_validation) <- colnames(prs_all_validation)
print(all.equal(test_validation,prs_all_validation))

Final_Coefficients_LDPred <- read.csv(paste0(trait,"_Final_Coefficients_LDPred2.csv"))
system(paste0("rm ",paste0(trait,"_Final_Coefficients_LDPred2.csv")))
prs.file <- Final_Coefficients_LDPred[,c("SNP","A1",paste0("LDPred2_SCORE",1:255,"_SUM"))]
write.table(prs.file,file = paste0(trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("plink2 --score ",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile all_chr --keep validation.txt --threads 1 --out test_validation_",trait))
test_validation <- read.delim(paste0("test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0(trait,"_prs_validation_ldpred2.sscore"), sep="")
system(paste0("rm ",paste0(trait,"_prs_validation_ldpred2.sscore")))
prs_all_validation <- prs_all_validation[,c(2,5:ncol(prs_all_validation))]
colnames(test_validation) <- colnames(prs_all_validation)
print(all.equal(test_validation,prs_all_validation))

Final_Coefficients_LASSOSUM <- read.csv(paste0(trait,"_Final_Coefficients_LASSOSum.csv"))
system(paste0("rm ",paste0(trait,"_Final_Coefficients_LASSOSum.csv")))
prs.file <- Final_Coefficients_LASSOSUM[,c("SNP","A1",paste0("LASSOSum2_SCORE",1:300,"_SUM"))]
write.table(prs.file,file = paste0(trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("plink2 --score ",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile all_chr --keep validation.txt --threads 1 --out test_validation_",trait))
test_validation <- read.delim(paste0("test_validation_",trait,".sscore"), header=FALSE, comment.char="#")
test_validation <- test_validation[,c(2,5:ncol(test_validation))]
prs_all_validation <- read.csv(paste0(trait,"_prs_validation_lassosum2.sscore"), sep="")
system(paste0("rm ",paste0(trait,"_prs_validation_lassosum2.sscore")))
prs_all_validation <- prs_all_validation[,c(2,5:ncol(prs_all_validation))]
colnames(test_validation) <- colnames(prs_all_validation)
print(all.equal(test_validation,prs_all_validation))



# ### Coefficient Way
# BMI_All <- cbind(BMI_Final_Coefficients_CT,BMI_Final_Coefficients_LDPred[,str_detect(colnames(BMI_Final_Coefficients_LDPred),"LDPred2")],BMI_Final_Coefficients_LASSOSUM[,str_detect(colnames(BMI_Final_Coefficients_LASSOSUM),"LASSOSum2")])
# 
# BMI_Final_Coefficients_SL <- read.csv("BMI_final_coef.csv")
# 
# CV_SL_Beta <- BMI_Final_Coefficients_CT[,c("CHR","SNP","REF","BP","A1")]
# CV_SL_Beta$BETA <- 0
# 
# for(i in 1:nrow(BMI_Final_Coefficients_SL)){
#   name_i <- BMI_Final_Coefficients_SL$Coef[i]
#   if(name_i %in% colnames(BMI_All)){
#     CV_SL_Beta$BETA <- CV_SL_Beta$BETA + BMI_All[,name_i]*BMI_Final_Coefficients_SL$Beta[i] 
#   }
# }
# 
# prs.file <- CV_SL_Beta[,c("SNP","A1","BETA")]
# write.table(prs.file,file = "BMI_Final_Score",col.names = T,row.names = F,quote=F)
# 
# system("plink2 --score BMI_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --bfile all_chr --keep validation.txt --threads 1 --out test_validation")
# 
# Best_Validation_All <- read.delim("BMI_Best_Validation_All.txt")
# test_validation <- read.delim("test_validation.sscore", header=FALSE, comment.char="#")
# 
# all.equal(Best_Validation_All$IID,test_validation[,2])
# all.equal(Best_Validation_All$prs,test_validation[,5])
# cor(Best_Validation_All$prs,test_validation[,5])


#### Linear Model Way
Final_Coefficients_SL <- read.csv(paste0(trait,"_Final_Coefficients_SL.csv"))
system(paste0("rm ",paste0(trait,"_Final_Coefficients_SL.csv")))

prs_tune_CT <- read.csv(paste0(trait,"_prs_all_tune.txt"), sep="")
prs_tune_LDPred2 <- read.delim(paste0(trait,"_prs_tune_ldpred2.sscore"))
prs_tune_LASSOSum2 <- read.delim(paste0(trait,"_prs_tune_lassosum2.sscore"))

system(paste0("rm ",paste0(trait,"_prs_all_tune.txt")))
system(paste0("rm ",paste0(trait,"_prs_tune_ldpred2.sscore")))
system(paste0("rm ",paste0(trait,"_prs_tune_lassosum2.sscore")))

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

## Merge covariates and y for tuning with the prs_mat
Best_Tune_All <- read.delim(paste0(trait,"_Best_Tune_All.txt"))
system(paste0("rm ",paste0(trait,"_Best_Tune_All.txt")))
Best_Tune_All <- left_join(Best_Tune_All,prs_tune_all,by = "IID")
Best_Tune_All <- Best_Tune_All[,colnames(Best_Tune_All) %in% c("prs",Final_Coefficients_SL$Coef)]

if(nrow(Final_Coefficients_SL) == 1){
  Beta_Star <- matrix(unname(coef(lm(prs ~.,data = Best_Tune_All))),ncol = 1)
  Beta_Star <- Beta_Star[-1,,drop = FALSE]
}else{
  Beta_Star <- matrix(unname(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All))),ncol = 1)
  Beta_Star <- Beta_Star[-1,,drop = FALSE] 
}

All <- cbind(Final_Coefficients_CT,Final_Coefficients_LDPred[,str_detect(colnames(Final_Coefficients_LDPred),"LDPred2")],Final_Coefficients_LASSOSUM[,str_detect(colnames(Final_Coefficients_LASSOSUM),"LASSOSum2")])
score_full <- All[,colnames(All) %in% names(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All)))[-1],drop = FALSE]
score_full <- score_full[,match(colnames(score_full),names(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All)))[-1]),drop = FALSE]
score_full <- as.matrix(score_full)

prs.file <- data.frame(SNP = All$SNP,A1 = All$A1,BETA = score_full %*% Beta_Star)
write.table(prs.file,file = paste0(trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("plink2 --score ",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --bfile all_chr --keep validation.txt --threads 1 --out test_validation_",trait))

system("rm all_chr.bed")
system("rm all_chr.bim")
system("rm all_chr.fam")
system("rm validation.txt")

Best_Validation_All <- read.delim(paste0(trait,"_Best_Validation_All.txt"))
system(paste0("rm ",paste0(trait,"_Best_Validation_All.txt")))
test_validation <- read.delim(paste0("test_validation_",trait,".sscore"), header=FALSE, comment.char="#")

all.equal(Best_Validation_All$IID,test_validation[,2])
all.equal(Best_Validation_All$prs,exp(matrix(unname(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All))),ncol = 1)[1] + test_validation[,5])/(1 + exp(matrix(unname(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All))),ncol = 1)[1] + test_validation[,5])))
cor(Best_Validation_All$prs,exp(matrix(unname(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All))),ncol = 1)[1] + test_validation[,5])/(1 + exp(matrix(unname(coef(lm(log(prs/(1 - prs)) ~.,data = Best_Tune_All))),ncol = 1)[1] + test_validation[,5])))

system(paste0("rm test_validation_",trait,".sscore"))
system(paste0("rm test_validation_",trait,".log"))

