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
rm(i)

write.table(prs_train_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_train_prs_all.txt",sep = "\t")
write.table(prs_tune_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_tune_prs_all.txt",sep = "\t")
write.table(prs_validation_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_validation_prs_all.txt",sep = "\t")

rm(list = ls())

## Train
for(i in 1:22){
  ## Train
  if(i == 1){
    prs_train_mat <- fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_train_chr",i,".sscore"))
  }else{
    prs_train_mat[,5:304] <- prs_train_mat[,5:304] + fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_train_chr",i,".sscore"))[,5:304]
  }
  
  ## Tune
  if(i == 1){
    prs_tune_mat <- fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_tune_chr",i,".sscore"))
  }else{
    prs_tune_mat[,5:304] <- prs_tune_mat[,5:304] + fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_tune_chr",i,".sscore"))[,5:304]
  }
  
  ## Validation
  if(i == 1){
    prs_validation_mat <- fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_validation_chr",i,".sscore"))
  }else{
    prs_validation_mat[,5:304] <- prs_validation_mat[,5:304] + fread(paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_validation_chr",i,".sscore"))[,5:304]
  }
}
rm(i)

write.table(prs_train_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_train_prs_all.txt",sep = "\t")
write.table(prs_tune_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_tune_prs_all.txt",sep = "\t")
write.table(prs_validation_mat,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_validation_prs_all.txt",sep = "\t")

# --------------------- Validation ---------------------
### LDPred2

rm(list=ls())
library(data.table)
library(dplyr)
library(bigsnpr)
library(boot)

h2_seq <- c(0.7, 1, 1.4) 
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

prs_mat_tune <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_tune_prs_all.txt")
prs_mat_validation <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_validation_prs_all.txt")

## Pull in Phenotypes/Covariates 
pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
colnames(pheno_tuning) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
## Merge covariates and y for tuning with the prs_mat
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,nrow(sets))
model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tuning)
for(k in 1:nrow(sets)){
  prs <- pheno_tuning[,paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad)
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
boot_r2 <- boot(data = data, statistic = R2Boot, R = 100000)

ci_result <- boot.ci(boot_r2, type = "bca")
r2.result <- data.frame(method = "LDPred2",
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

prs_mat_tune <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_tune_prs_all.txt")
prs_mat_validation <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_validation_prs_all.txt")

## Pull in Phenotypes/Covariates 
pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
colnames(pheno_tuning) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
## Merge covariates and y for tuning with the prs_mat
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
colnames(pheno_vad) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,300)
model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tuning)
for(k in 1:300){
  prs <- pheno_tuning[,paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

model.vad.null  <-  lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_vad)
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
boot_r2 <- boot(data = data, statistic = R2Boot, R = 100000)

ci_result <- boot.ci(boot_r2, type = "bca")
r2.result <- data.frame(method = "LASSOSUM2",
                        r2 = r2,
                        r2_low = ci_result$bca[4],
                        r2_high = ci_result$bca[5]
)

## Save the R2 for the validation set w/ its confidence bounds, as well as the R2 tuning vector
lassosum2.result <- list(r2.result,r2_tun_vec)
rm(list=setdiff(ls(), c("lassosum2.result","ldpred2.result")))

save(ldpred2.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_result.RData")
save(lassosum2.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2_result.RData")