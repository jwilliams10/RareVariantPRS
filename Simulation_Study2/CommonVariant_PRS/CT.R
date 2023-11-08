rm(list = ls())
library(data.table)
library(dplyr)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

dat <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/GWAS_Summary_Statistics/Y_Train",i,".Y.glm.linear"), header=FALSE, comment.char="#")
colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
dat <- dat[dat$TEST == "ADD",]

dat <- dat[,c("CHROM","ID","REF","POS","A1","BETA","P")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1","BETA","P")

write.table(dat,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/GWAS_Summary_Statistics/tmp",i,".txt"),col.names = T,row.names = F,quote=F)

pthr <- 1
r2thr <- 0.1
kbpthr <- 500
# Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/reference_CT.txt --clump /data/williamsjacr/UKB_WES_Simulation/Simulation2/GWAS_Summary_Statistics/tmp",i,".txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ","/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/LD_clump",i))

################################################

## This is the beginning of the thresholding step of clumping and threshold

#p-value thresholds
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
LD <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/LD_clump",i,".clumped")))
# grab the index SNP for each clump
clump.snp <- LD[,3,drop=F]
# join against the summary data, this is now a dataset with n = number of index SNPs
prs.all <- left_join(clump.snp,dat)

n_pthres <- length(pthres)

## Write Coefficients of index SNPs to use later
prs.file <- prs.all[,c("SNP","A1","BETA")]
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_coeff",i),col.names = T,row.names = F,quote=F)

# Write p-values to file to use later
p.value.file <- prs.all[,c("SNP","P")]
write.table(p.value.file,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/p_value",i),col.names = T,row.names = F,quote=F)

## Write a file with the p-value thresholds 
q_range <- data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
temp <- 1
for(k in 1:length(pthres)){
  q_range[temp,1] <- paste0("p_value_",k)
  q_range[temp,3] <- pthres[k]
  temp <- temp + 1
}
q_range <- q_range[1:(temp-1),]
write.table(q_range,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/q_range_file",i),row.names = F,col.names = F,quote=F)

#PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
#PRS = G*beta
# Literally just multiply the score file (weighted or unweighted coefficients) by the G matrix, q-score-range is only for C + T, for LD pred score file would be weight coefficients 
system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/q_range_file",i," /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/p_value",i," header --threads 2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_coeff",i," cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/train.txt --out /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_train",i))
system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/q_range_file",i," /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/p_value",i," header --threads 2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_coeff",i," cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/tune.txt --out /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_tune",i))
system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/q_range_file",i," /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/p_value",i," header --threads 2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_coeff",i," cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/validation.txt --out /data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_validation",i))

#########################################################################

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_train",i,".p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_train <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_train)[2] <- "id"

write.table(prs_mat_train,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_train",i,".txt"),row.names = F)

prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_tune",i,".p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_tune <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_tune)[2] <- "id"

write.table(prs_mat_tune,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_tune",i,".txt"),row.names = F)

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_validation",i,".p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_validation <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_validation)[2] <- "id"

write.table(prs_mat_validation,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_all_validation",i,".txt"),row.names = F)


## Pull in Phenotypes/Covariates 

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("id","Y")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "id")

pheno_tuning <- Y_tune[[i]]
colnames(pheno_tuning) <- c("id","Y")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "id")

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("id","Y")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "id")

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
#write down best prs in the prs folder/ save it and store it
prs_train_max <- pheno_train[,c("id",paste0("p_value_",idx))]
colnames(prs_train_max) <- c("IID","prs")
write.table(prs_train_max, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_train_best",i,".txt"),row.names = F)

prs_tune_max <- pheno_tuning[,c("id",paste0("p_value_",idx))]
colnames(prs_tune_max) <- c("IID","prs")
write.table(prs_tune_max, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_tune_best",i,".txt"),row.names = F)

prs_vad_max <- pheno_vad[,c("id",paste0("p_value_",idx))]
colnames(prs_vad_max) <- c("IID","prs")
write.table(prs_vad_max, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/prs_validation_best",i,".txt"),row.names = F)

#evaluate the best threshold based on the tuning on the validation dataset
model.vad.null  <-  lm(Y~1,data=pheno_vad)
prs <- pheno_vad[,paste0("p_value_",idx)]
model.vad.prs <- lm(model.vad.null$residual~prs)
r2 <- summary(model.vad.prs)$r.square

## bootstrap the R2 to provide an approximate distribution 
data <- data.frame(y = model.vad.null$residual, x = prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
library(boot)
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
r2.result <- data.frame(method = "CT",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
ct.result <- list(r2.result,r2_tun_vec)
save(ct.result, file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result",i,".RData"))

