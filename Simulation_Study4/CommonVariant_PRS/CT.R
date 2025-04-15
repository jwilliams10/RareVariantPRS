rm(list = ls())
library(data.table)
library(dplyr)
library(boot)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

dat <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/GWAS_Summary_Statistics/Y_Train",i,".Y.glm.linear"), header=FALSE, comment.char="#")
colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
dat <- dat[dat$TEST == "ADD",]

dat <- dat[,c("CHROM","ID","REF","POS","A1","BETA","P")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1","BETA","P")

write.table(dat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/GWAS_Summary_Statistics/tmp",i,".txt"),col.names = T,row.names = F,quote=F)

pthr <- 1
r2thr <- 0.1
kbpthr <- 500
# Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/reference_CT.txt --clump /data/williamsjacr/UKB_WES_Simulation/Simulation4/GWAS_Summary_Statistics/tmp",i,".txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ","/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/LD_clump",i))

system(paste0("rm /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/LD_clump",i,".log"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/GWAS_Summary_Statistics/tmp",i,".txt")))

################################################

## This is the beginning of the thresholding step of clumping and threshold

#p-value thresholds
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
LD <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/LD_clump",i,".clumped")))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/LD_clump",i,".clumped")))
# grab the index SNP for each clump
clump.snp <- LD[,3,drop=F]
# join against the summary data, this is now a dataset with n = number of index SNPs
prs.all <- left_join(clump.snp,dat)

n_pthres <- length(pthres)

## Write Coefficients of index SNPs to use later
prs.file <- prs.all[,c("SNP","A1","BETA")]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_coeff",i),col.names = T,row.names = F,quote=F)

# Write p-values to file to use later
p.value.file <- prs.all[,c("SNP","P")]
write.table(p.value.file,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/p_value",i),col.names = T,row.names = F,quote=F)

## Write a file with the p-value thresholds 
q_range <- data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
temp <- 1
for(k in 1:length(pthres)){
  q_range[temp,1] <- paste0("p_value_",k)
  q_range[temp,3] <- pthres[k]
  temp <- temp + 1
}
q_range <- q_range[1:(temp-1),]
write.table(q_range,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/q_range_file",i),row.names = F,col.names = F,quote=F)

#PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
#PRS = G*beta
# Literally just multiply the score file (weighted or unweighted coefficients) by the G matrix, q-score-range is only for C + T, for LD pred score file would be weight coefficients 
system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/q_range_file",i," /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/p_value",i," header --threads 2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_coeff",i," cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/train.txt --out /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_train",i))
system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/q_range_file",i," /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/p_value",i," header --threads 2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_coeff",i," cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/tune.txt --out /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_tune",i))
system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/q_range_file",i," /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/p_value",i," header --threads 2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_coeff",i," cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/validation.txt --out /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_validation",i))

system(paste0("rm /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_train",i,".log"))
system(paste0("rm /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_tune",i,".log"))
system(paste0("rm /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_validation",i,".log"))
system(paste0("rm /data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/q_range_file",i))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/p_value",i)))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_coeff",i)))

#########################################################################

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_train",i,".p_value_",k,".sscore"))
  system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_train",i,".p_value_",k,".sscore")))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_train <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_train)[2] <- "IID"

write.table(prs_mat_train,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_all_train",i,".txt"),row.names = F)

prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_tune",i,".p_value_",k,".sscore"))
  system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_tune",i,".p_value_",k,".sscore")))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_tune <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_tune)[2] <- "IID"

write.table(prs_mat_tune,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_all_tune",i,".txt"),row.names = F)

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_validation",i,".p_value_",k,".sscore"))
  system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_validation",i,".p_value_",k,".sscore")))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_validation <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_validation)[2] <- "IID"

write.table(prs_mat_validation,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/prs_all_validation",i,".txt"),row.names = F)


## Pull in Phenotypes/Covariates 

load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- Y_tune[[i]]
colnames(pheno_tuning) <- c("IID","Y")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,length(pthres))
model.null <- lm(Y~1,data=pheno_tuning)
for(k in 1:length(pthres)){
  prs <- pheno_tuning[,paste0("p_value_",k)]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}
#find best p-value threshold
idx <- which.max(r2_tun_vec)





load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad$y_validation <- NA
pheno_vad$y_validation <- lm(Y~1,data=pheno_vad)$residual

pheno_validation_raw <- pheno_vad
pheno_validation_adjusted <- pheno_vad
pheno_validation_adjusted <- inner_join(pheno_validation_adjusted,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])

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

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]

pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
pheno_validation_raw_AMR$y_validation <- scale(pheno_validation_raw_AMR$y_validation)
pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)

pheno_validation_raw_EUR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EUR[,paste0("p_value_",idx)])
pheno_validation_raw_SAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_SAS[,paste0("p_value_",idx)])
pheno_validation_raw_AMR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AMR[,paste0("p_value_",idx)])
pheno_validation_raw_AFR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AFR[,paste0("p_value_",idx)])

Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = boot_data))[2]
  return(c(result))
}

R2_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = boot_data))$r.squared
  return(c(result))
}

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_EUR))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EUR <- sd(boot_R2$t)
R2_lower_validation_raw_EUR <- R2_ci$basic[4]
R2_upper_validation_raw_EUR <- R2_ci$basic[5]

beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_SAS))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_SAS <- sd(boot_R2$t)
R2_lower_validation_raw_SAS <- R2_ci$basic[4]
R2_upper_validation_raw_SAS <- R2_ci$basic[5]

beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AMR))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AMR <- sd(boot_R2$t)
R2_lower_validation_raw_AMR <- R2_ci$basic[4]
R2_upper_validation_raw_AMR <- R2_ci$basic[5]

beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AFR))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AFR <- sd(boot_R2$t)
R2_lower_validation_raw_AFR <- R2_ci$basic[4]
R2_upper_validation_raw_AFR <- R2_ci$basic[5]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_EUR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
R2_lower_validation_adjusted_EUR <- R2_ci$basic[4]
R2_upper_validation_adjusted_EUR <- R2_ci$basic[5]

beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_SAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_SAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_SAS <- R2_ci$basic[5]

beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AMR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AMR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AMR <- R2_ci$basic[5]

beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AFR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AFR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AFR <- R2_ci$basic[5]

CT_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR"), 
                         beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR), 
                         beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR), 
                         beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR), 
                         beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR), 
                         R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR),
                         R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR),
                         R2_lower_raw = c(R2_lower_validation_raw_EUR,R2_lower_validation_raw_SAS,R2_lower_validation_raw_AMR,R2_lower_validation_raw_AFR),
                         R2_upper_raw = c(R2_upper_validation_raw_EUR,R2_upper_validation_raw_SAS,R2_upper_validation_raw_AMR,R2_upper_validation_raw_AFR),
                         beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR), 
                         beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR), 
                         beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR), 
                         beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR), 
                         R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR),
                         R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR),
                         R2_lower_adjusted = c(R2_lower_validation_adjusted_EUR,R2_lower_validation_adjusted_SAS,R2_lower_validation_adjusted_AMR,R2_lower_validation_adjusted_AFR),
                         R2_upper_adjusted = c(R2_upper_validation_adjusted_EUR,R2_upper_validation_adjusted_SAS,R2_upper_validation_adjusted_AMR,R2_upper_validation_adjusted_AFR))

write.csv(CT_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/Best_Betas",i,".csv"),row.names = FALSE)