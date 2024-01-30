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
  system("rm regenie_step2_act_Asthma.regenie")
}else if(trait == 2){
  trait <- "CAD"
  dat <- read.csv("regenie_step2_act_CAD.regenie", sep="")
  system("rm regenie_step2_act_CAD.regenie")
}else if(trait == 3){
  trait <- "T2D"
  dat <- read.csv("regenie_step2_act_T2D.regenie", sep="")
  system("rm regenie_step2_act_T2D.regenie")
}else if(trait == 4){
  trait <- "Breast"
  dat <- read.csv("regenie_step2_bp_Breast.regenie", sep="")
  system("rm regenie_step2_bp_Breast.regenie")
}else{
  trait <- "Prostate"
  dat <- read.csv("regenie_step2_bp_Prostate.regenie", sep="")
  system("rm regenie_step2_bp_Prostate.regenie")
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

load("all_phenotypes.RData")

system(paste0("rm All_Train.txt"))
system(paste0("rm All_Tune.txt"))
system(paste0("rm All_Validation.txt"))
system(paste0("rm all_phenotypes.RData"))

pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

#calculate AUC for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
AUC_tun_vec <- rep(0,length(pthres))
if(trait %in% c("Breast","Prostate")){
  confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
}else{
  confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
}

for(k in 1:length(pthres)){
  d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",k))]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",k),
                        confounders = confounders,
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

#evaluate the best threshold based on the tuning on the validation dataset
d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("p_value_",idx),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "CT_EUR",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ct.result <- list(AUC.result,AUC_tun_vec)
save(ct.result, file = paste0(trait,"_CT_result_EUR.RData")) 


#evaluate the best threshold based on the tuning on the validation dataset
d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("p_value_",idx),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "CT_NonEur",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ct.result <- list(AUC.result,AUC_tun_vec)
save(ct.result, file = paste0(trait,"_CT_result_NonEur.RData")) 


#evaluate the best threshold based on the tuning on the validation dataset
d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("p_value_",idx),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "CT_UNK",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ct.result <- list(AUC.result,AUC_tun_vec)
save(ct.result, file = paste0(trait,"_CT_result_UNK.RData")) 


#evaluate the best threshold based on the tuning on the validation dataset
d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("p_value_",idx),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "CT_SAS",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ct.result <- list(AUC.result,AUC_tun_vec)
save(ct.result, file = paste0(trait,"_CT_result_SAS.RData")) 


#evaluate the best threshold based on the tuning on the validation dataset
d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("p_value_",idx),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "CT_MIX",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ct.result <- list(AUC.result,AUC_tun_vec)
save(ct.result, file = paste0(trait,"_CT_result_MIX.RData"))


#evaluate the best threshold based on the tuning on the validation dataset
d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("p_value_",idx),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "CT_AFR",
                         AUC = AUC,
                         AUC_low = ci_result$percent[4],
                         AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ct.result <- list(AUC.result,AUC_tun_vec)
save(ct.result, file = paste0(trait,"_CT_result_AFR.RData"))


#evaluate the best threshold based on the tuning on the validation dataset
if(trait %in% c("Prostate","CAD")){
  ct.result <- NA
}else{
  d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("p_value_",idx),
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("p_value_",idx),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  
  AUC.result <- data.frame(method = "CT_EAS",
                           AUC = AUC,
                           AUC_low = ci_result$percent[4],
                           AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result <- list(AUC.result,AUC_tun_vec) 
}
save(ct.result, file = paste0(trait,"_CT_result_EAS.RData"))
