rm(list = ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)

## Train
for(i in 1:22){
  ## Train
  if(i == 1){
    prs_train_mat <- fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_train_chr",i,".sscore"))
  }else{
    prs_train_mat[,5:55] <- prs_train_mat[,5:55] + fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_train_chr",i,".sscore"))[,5:55]
  }
  
  ## Tune
  if(i == 1){
    prs_tune_mat <- fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_tune_chr",i,".sscore"))
  }else{
    prs_tune_mat[,5:55] <- prs_tune_mat[,5:55] + fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_tune_chr",i,".sscore"))[,5:55]
  }
  
  ## Validation
  if(i == 1){
    prs_validation_mat <- fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_validation_chr",i,".sscore"))
  }else{
    prs_validation_mat[,5:55] <- prs_validation_mat[,5:55] + fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_validation_chr",i,".sscore"))[,5:55]
  }
}
# no NA's
# apply(prs_mat,2,function(x){sum(is.na(x))})
rm(i)

write.table(prs_train_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_train_prs_all.txt",sep = "\t")
write.table(prs_tune_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_tune_prs_all.txt",sep = "\t")
write.table(prs_validation_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_validation_prs_all.txt",sep = "\t")

rm(list = ls())
for(i in 1:22){
  if(i == 1){
    prs_mat <- fread(paste0("/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/lassosum2/score/lassosum2-chr",i,".sscore"))
  }else{
    prs_mat[,5:304] <- prs_mat[,5:304] + fread(paste0("/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/lassosum2/score/lassosum2-chr",i,".sscore"))[,5:304]
  }
}
# no NA's
# apply(prs_mat,2,function(x){sum(is.na(x))})
rm(i)

write.table(prs_mat,file="/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/lassosum2/lassosum2_prs_all.txt",sep = "\t")

# --------------------- Validation ---------------------
### LDPred2

rm(list=ls())
library(data.table)
library(dplyr)
library(bigsnpr)
library(boot)
race = c('EUR','AFR','AMR')[2]
trait = c("HDL","LDL","logTG","TC")[1]

h2_seq <- c(0.7, 1, 1.4) 
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

prs_mat <- read.delim("/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/ldpred2/ldpred2_prs_all.txt")

## Pull in Phenotypes/Covariates 
pheno.dir <- "/data/BB_Bioinformatics/ProjectData/UKBB_sub/phenotype/"
pheno_tuning <- as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",race,"_tuning.txt")))
pheno_tuning <- pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",race,"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) <- c('IID','y','sex','age',paste0('pc',1:10))

## So there was 10,000 original but only 8,742 with complete cases
pheno_tuning_com <- pheno_tuning[complete.cases(pheno_tuning$y),]

## Merge covariates and y for tuning with the prs_mat
pheno_tuning <- left_join(pheno_tuning_com,prs_mat,by = "IID")

pheno_vad <- as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",race,"_validation.txt")))
pheno_vad <- pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) <- c('IID','y','sex','age',paste0('pc',1:10))

## Was 10,000 then dropped to 8,715 with complete cases
pheno_vad_com <- pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad <- left_join(pheno_vad_com,prs_mat,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,nrow(sets))
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
for(k in 1:nrow(sets)){
  prs <- pheno_tuning[,paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

model.vad.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
prs <- pheno_vad[,paste0("SCORE",idx,"_SUM")]
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
boot_r2 <- boot(data = data, statistic = R2Boot, R = 10000)

ci_result <- boot.ci(boot_r2, type = "bca")
r2.result <- data.frame(eth = race,
                        trait = trait,
                        method = "CT",
                        r2 = r2,
                        r2_low = ci_result$bca[4],
                        r2_high = ci_result$bca[5]
)

## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
ldpred2.result <- list(r2.result,r2_tun_vec)

##################################################################################

rm(list=setdiff(ls(), "ldpred2.result"))
library(data.table)
library(dplyr)
library(bigsnpr)
library(boot)
race = c('EUR','AFR','AMR')[2]
trait = c("HDL","LDL","logTG","TC")[1]

prs_mat <- read.delim("/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/lassosum2/lassosum2_prs_all.txt")

## Pull in Phenotypes/Covariates 
pheno.dir <- "/data/BB_Bioinformatics/ProjectData/UKBB_sub/phenotype/"
pheno_tuning <- as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",race,"_tuning.txt")))
pheno_tuning <- pheno_tuning[,1:2]
covar <- as.data.frame(fread(paste0(pheno.dir,"/covariates/tuning+validation/",race,"_all_data.txt")))
pheno_tuning <- left_join(pheno_tuning, covar)
colnames(pheno_tuning) <- c('IID','y','sex','age',paste0('pc',1:10))

## So there was 10,000 original but only 8,742 with complete cases
pheno_tuning_com <- pheno_tuning[complete.cases(pheno_tuning$y),]

## Merge covariates and y for tuning with the prs_mat
pheno_tuning <- left_join(pheno_tuning_com,prs_mat,by = "IID")

pheno_vad <- as.data.frame(fread(paste0(pheno.dir,trait,"/tuning+validation/",race,"_validation.txt")))
pheno_vad <- pheno_vad[,1:2]
pheno_vad <- left_join(pheno_vad, covar)
colnames(pheno_vad) <- c('IID','y','sex','age',paste0('pc',1:10))

## Was 10,000 then dropped to 8,715 with complete cases
pheno_vad_com <- pheno_vad[complete.cases(pheno_vad$y),]
pheno_vad <- left_join(pheno_vad_com,prs_mat,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,300)
model.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_tuning)
for(k in 1:300){
  prs <- pheno_tuning[,paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

model.vad.null <- lm(y~pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+age+sex,data=pheno_vad)
prs <- pheno_vad[,paste0("SCORE",idx,"_SUM")]
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
boot_r2 <- boot(data = data, statistic = R2Boot, R = 10000)

ci_result <- boot.ci(boot_r2, type = "bca")
r2.result <- data.frame(eth = race,
                        trait = trait,
                        method = "CT",
                        r2 = r2,
                        r2_low = ci_result$bca[4],
                        r2_high = ci_result$bca[5]
)

## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
lassosum2.result <- list(r2.result,r2_tun_vec)
rm(list=setdiff(ls(), c("lassosum2.result","ldpred2.result")))