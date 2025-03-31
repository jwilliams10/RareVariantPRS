rm(list = ls())

if(!("data.table" %in% rownames(installed.packages()))){
  install.packages("data.table",quiet = TRUE)
}

if(!("dplyr" %in% rownames(installed.packages()))){
  install.packages("dplyr",quiet = TRUE)
}

if(!("boot" %in% rownames(installed.packages()))){
  install.packages("boot",quiet = TRUE)
}
if(!("RISCA" %in% rownames(installed.packages()))){
  install.packages("RISCA",quiet = TRUE)
}

library(data.table)
library(dplyr)
library(RISCA)
library(boot)

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/CommonVariant_PRS/CT_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/CommonVariant_PRS/CT_Binary.sh -icmd="bash CT_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/CT/ --priority low --instance-type mem1_ssd1_v2_x36
# done

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "Asthma"
  dat <- read.csv("regenie_step2_act_Asthma.regenie", sep="")
}else if(trait == 2){
  trait <- "CAD"
  dat <- read.csv("regenie_step2_act_CAD.regenie", sep="")
}else if(trait == 3){
  trait <- "T2D"
  dat <- read.csv("regenie_step2_act_T2D.regenie", sep="")
}else if(trait == 4){
  trait <- "Breast"
  dat <- read.csv("regenie_step2_bp_Breast.regenie", sep="")
}else{
  trait <- "Prostate"
  dat <- read.csv("regenie_step2_bp_Prostate.regenie", sep="")
}

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
dat$P <- 10^(-1*dat$LOG10P)

dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1","BETA","P")

write.table(dat,file = paste0("",trait,"_assoc.txt"),col.names = T,row.names = F,quote=F)

pthr <- 1
r2thr <- 0.1
kbpthr <- 500
# Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
system(paste0("plink --bfile all_chr_reference --clump ",trait,"_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",trait,"_LDclump"))

system(paste0("rm all_chr_reference.bed"))
system(paste0("rm all_chr_reference.bim"))
system(paste0("rm all_chr_reference.fam"))
system(paste0("rm ",trait,"_assoc.txt"))

system(paste0("rm ",trait,"_LDclump.log"))
system(paste0("rm ",trait,"_LDclump.nosex"))

################################################

## This is the beginning of the thresholding step of clumping and threshold

#p-value thresholds
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
LD <- as.data.frame(fread(paste0(trait,"_LDclump.clumped")))
system(paste0("rm ",trait,"_LDclump.clumped"))
# grab the index SNP for each clump
clump.snp <- LD[,3,drop=F]
# join against the summary data, this is now a dataset with n = number of index SNPs
prs.all <- left_join(clump.snp,dat)

n_pthres <- length(pthres)

## Write Coefficients of index SNPs to use later
prs.file <- prs.all[,c("SNP","A1","BETA")]
write.table(prs.file,file = paste0(trait,"_prs_coeff"),col.names = T,row.names = F,quote=F)

# Write p-values to file to use later
p.value.file <- prs.all[,c("SNP","P")]
write.table(p.value.file,file = paste0(trait,"_p_value"),col.names = T,row.names = F,quote=F)

beta_final <- inner_join(prs.file,p.value.file)

## Write a file with the p-value thresholds 
q_range <- data.frame(a = paste0("p_value_",1:length(pthres)),b = 0,c = pthres)

write.table(q_range,file = "q_range_file",row.names = F,col.names = F,quote=F)

#PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
#PRS = G*beta
# Literally just multiply the score file (weighted or unweighted coefficients) by the G matrix, q-score-range is only for C + T, for LD pred score file would be weight coefficients 

system(paste0("plink2 --q-score-range q_range_file ",trait,"_p_value header --threads 2 --score ",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile all_chr --keep train.txt --out ",trait,"_prs_train"))
system(paste0("plink2 --q-score-range q_range_file ",trait,"_p_value header --threads 2 --score ",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile all_chr --keep tune.txt --out ",trait,"_prs_tune"))
system(paste0("plink2 --q-score-range q_range_file ",trait,"_p_value header --threads 2 --score ",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile all_chr --keep validation.txt --out ",trait,"_prs_validation"))
#########################################################################

system(paste0("rm ",trait,"_prs_validation.log"))
system(paste0("rm ",trait,"_prs_tune.log"))
system(paste0("rm ",trait,"_prs_train.log"))

system(paste0("rm q_range_file"))
system(paste0("rm ",trait,"_p_value"))
system(paste0("rm ",trait,"_prs_coeff"))

system(paste0("rm all_chr.bed"))
system(paste0("rm all_chr.bim"))
system(paste0("rm all_chr.fam"))

system(paste0("rm train.txt"))
system(paste0("rm tune.txt"))
system(paste0("rm validation.txt"))

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0(trait,"_prs_train.p_value_",k,".sscore"))
  system(paste0("rm ",trait,"_prs_train.p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_train <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_train)[2] <- "IID"

write.table(prs_mat_train,file = paste0(trait,"_prs_all_train.txt"),row.names = F)

prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0(trait,"_prs_tune.p_value_",k,".sscore"))
  system(paste0("rm ",trait,"_prs_tune.p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_tune <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_tune)[2] <- "IID"

write.table(prs_mat_tune,file = paste0(trait,"_prs_all_tune.txt"),row.names = F)

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0(trait,"_prs_validation.p_value_",k,".sscore"))
  system(paste0("rm ",trait,"_prs_validation.p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_validation <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_validation)[2] <- "IID"

write.table(prs_mat_validation,file = paste0(trait,"_prs_all_validation.txt"),row.names = F)


## Pull in Phenotypes/Covariates 
pheno_train <- read.delim("All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate AUC for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
AUC_tun_vec <- rep(0,length(pthres))
if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <- paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}

for(k in 1:length(pthres)){
  d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",k))]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",k),
                        confounders = as.formula(confounders),
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  
  AUC_tun_vec[k] <- roc_obj$auc
}

#find best p-value threshold
idx <- which.max(AUC_tun_vec)
#write down best prs in the prs folder/ save it and store it
prs_train_max <- pheno_train[,c("IID","FID",paste0("p_value_",idx))]
colnames(prs_train_max) <- c("IID","FID","prs")
write.table(prs_train_max, file = paste0(trait,"_prs_train_best.txt"),row.names = F)

prs_tune_max <- pheno_tuning[,c("IID","FID",paste0("p_value_",idx))]
colnames(prs_tune_max) <- c("IID","FID","prs")
write.table(prs_tune_max, file = paste0(trait,"_prs_tune_best.txt"),row.names = F)

prs_vad_max <- pheno_vad[,c("IID","FID",paste0("p_value_",idx))]
colnames(prs_vad_max) <- c("IID","FID","prs")
write.table(prs_vad_max, file = paste0(trait,"_prs_validation_best.txt"),row.names = F)

##### Final Coefficients
if(trait == "Asthma"){
  dat <- read.csv("regenie_step2_act_Asthma.regenie", sep="")
  system("rm regenie_step2_act_Asthma.regenie")
}else if(trait == "CAD"){
  dat <- read.csv("regenie_step2_act_CAD.regenie", sep="")
  system("rm regenie_step2_act_CAD.regenie")
}else if(trait == "T2D"){
  dat <- read.csv("regenie_step2_act_T2D.regenie", sep="")
  system("rm regenie_step2_act_T2D.regenie")
}else if(trait == "Breast"){
  dat <- read.csv("regenie_step2_bp_Breast.regenie", sep="")
  system("rm regenie_step2_bp_Breast.regenie")
}else{
  dat <- read.csv("regenie_step2_bp_Prostate.regenie", sep="")
  system("rm regenie_step2_bp_Prostate.regenie")
}

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
dat$P <- 10^(-1*dat$LOG10P)

dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1","BETA","P")

dat$BETA[is.na(dat$P)] <- 0
dat$P[is.na(dat$P)] <- 0

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

dat <- dat[,c("CHR","SNP","REF","BP","A1","P",paste0("CT_p_value_",1:9))]

write.csv(dat,file = paste0(trait,"_Final_Coefficients.csv"),row.names = FALSE)





load("all_phenotypes.RData")

system(paste0("rm All_Train.txt"))
system(paste0("rm All_Tune.txt"))
system(paste0("rm All_Validation.txt"))
system(paste0("rm all_phenotypes.RData"))

pheno_validation_raw <- pheno_vad
pheno_validation_adjusted <- pheno_vad

mod <- lm(as.formula(paste0(paste0("p_value_",idx),"~pc1 + pc2 + pc3 + pc4 + pc5")),data = pheno_validation_adjusted)
R <- mod$residuals
tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
y_hat <- predict(mod,tmp)
if(sum(sqrt(y_hat)) == 0){
  pheno_validation_adjusted[,paste0("p_value_",idx)] <- 0
}else{
  pheno_validation_adjusted[,paste0("p_value_",idx)] <- R/sqrt(y_hat)
}

pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EUR[,paste0("p_value_",idx)])
pheno_validation_raw_SAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_SAS[,paste0("p_value_",idx)])
pheno_validation_raw_AMR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AMR[,paste0("p_value_",idx)])
pheno_validation_raw_AFR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AFR[,paste0("p_value_",idx)])
pheno_validation_raw_EAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EAS[,paste0("p_value_",idx)])

Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
  return(c(result))
}

AUC_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  return(c(result))
}

beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

auc_validation_raw_EUR <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_raw_EUR[!is.na(pheno_validation_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_EUR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_EUR <- sd(boot_auc$t)
auc_lower_validation_raw_EUR <- auc_ci$basic[4]
auc_upper_validation_raw_EUR <- auc_ci$basic[5]

beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

auc_validation_raw_SAS <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_raw_SAS[!is.na(pheno_validation_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_SAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_SAS <- sd(boot_auc$t)
auc_lower_validation_raw_SAS <- auc_ci$basic[4]
auc_upper_validation_raw_SAS <- auc_ci$basic[5]

beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

auc_validation_raw_AMR <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_raw_AMR[!is.na(pheno_validation_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_AMR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_AMR <- sd(boot_auc$t)
auc_lower_validation_raw_AMR <- auc_ci$basic[4]
auc_upper_validation_raw_AMR <- auc_ci$basic[5]

beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

auc_validation_raw_AFR <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_raw_AFR[!is.na(pheno_validation_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_AFR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_AFR <- sd(boot_auc$t)
auc_lower_validation_raw_AFR <- auc_ci$basic[4]
auc_upper_validation_raw_AFR <- auc_ci$basic[5]

beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EAS <- sd(boot_beta$t)
beta_lower_validation_raw_EAS <- beta_ci$basic[4]
beta_upper_validation_raw_EAS <- beta_ci$basic[5]

auc_validation_raw_EAS <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_raw_EAS[!is.na(pheno_validation_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_raw_EAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_raw_EAS <- sd(boot_auc$t)
auc_lower_validation_raw_EAS <- auc_ci$basic[4]
auc_upper_validation_raw_EAS <- auc_ci$basic[5]

beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

auc_validation_adjusted_EUR <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_adjusted_EUR[!is.na(pheno_validation_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_EUR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_EUR <- sd(boot_auc$t)
auc_lower_validation_adjusted_EUR <- auc_ci$basic[4]
auc_upper_validation_adjusted_EUR <- auc_ci$basic[5]

beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

auc_validation_adjusted_SAS <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_adjusted_SAS[!is.na(pheno_validation_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_SAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_SAS <- sd(boot_auc$t)
auc_lower_validation_adjusted_SAS <- auc_ci$basic[4]
auc_upper_validation_adjusted_SAS <- auc_ci$basic[5]

beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

auc_validation_adjusted_AMR <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_adjusted_AMR[!is.na(pheno_validation_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_AMR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_AMR <- sd(boot_auc$t)
auc_lower_validation_adjusted_AMR <- auc_ci$basic[4]
auc_upper_validation_adjusted_AMR <- auc_ci$basic[5]

beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

auc_validation_adjusted_AFR <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_adjusted_AFR[!is.na(pheno_validation_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_AFR, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_AFR <- sd(boot_auc$t)
auc_lower_validation_adjusted_AFR <- auc_ci$basic[4]
auc_upper_validation_adjusted_AFR <- auc_ci$basic[5]

beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~",paste0("p_value_",idx),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]

auc_validation_adjusted_EAS <- roc.binary(status = trait,variable = paste0("p_value_",idx),confounders = as.formula(confounders),data = pheno_validation_adjusted_EAS[!is.na(pheno_validation_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_auc <- boot(data = pheno_validation_adjusted_EAS, statistic = AUC_Boot, R = 1000)
auc_ci <- boot.ci(boot_auc, type = "basic")
auc_se_validation_adjusted_EAS <- sd(boot_auc$t)
auc_lower_validation_adjusted_EAS <- auc_ci$basic[4]
auc_upper_validation_adjusted_EAS <- auc_ci$basic[5]

CT_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                         beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                         beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR,beta_lower_validation_raw_EAS), 
                         beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR,beta_upper_validation_raw_EAS), 
                         AUC_raw = c(auc_validation_raw_EUR,auc_validation_raw_SAS,auc_validation_raw_AMR,auc_validation_raw_AFR,auc_validation_raw_EAS),
                         AUC_se_raw = c(auc_se_validation_raw_EUR,auc_se_validation_raw_SAS,auc_se_validation_raw_AMR,auc_se_validation_raw_AFR,auc_se_validation_raw_EAS),
                         AUC_lower_raw = c(auc_lower_validation_raw_EUR,auc_lower_validation_raw_SAS,auc_lower_validation_raw_AMR,auc_lower_validation_raw_AFR,auc_lower_validation_raw_EAS),
                         AUC_upper_raw = c(auc_upper_validation_raw_EUR,auc_upper_validation_raw_SAS,auc_upper_validation_raw_AMR,auc_upper_validation_raw_AFR,auc_upper_validation_raw_EAS),
                         beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                         beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                         beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR,beta_lower_validation_adjusted_EAS), 
                         beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR,beta_upper_validation_adjusted_EAS), 
                         AUC_adjusted = c(auc_validation_adjusted_EUR,auc_validation_adjusted_SAS,auc_validation_adjusted_AMR,auc_validation_adjusted_AFR,auc_validation_adjusted_EAS),
                         AUC_se_adjusted = c(auc_se_validation_adjusted_EUR,auc_se_validation_adjusted_SAS,auc_se_validation_adjusted_AMR,auc_se_validation_adjusted_AFR,auc_se_validation_adjusted_EAS),
                         AUC_lower_adjusted = c(auc_lower_validation_adjusted_EUR,auc_lower_validation_adjusted_SAS,auc_lower_validation_adjusted_AMR,auc_lower_validation_adjusted_AFR,auc_lower_validation_adjusted_EAS),
                         AUC_upper_adjusted = c(auc_upper_validation_adjusted_EUR,auc_upper_validation_adjusted_SAS,auc_upper_validation_adjusted_AMR,auc_upper_validation_adjusted_AFR,auc_upper_validation_adjusted_EAS))

write.csv(CT_Results,file = paste0(trait,"Best_Betas.csv"),row.names = FALSE)