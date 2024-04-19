rm(list=ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)
library(boot)
library(bigstatsr)
library(RISCA)

snp_subset <- function(x, ind.row = bigstatsr::rows_along(x$genotypes), ind.col = bigstatsr::cols_along(x$genotypes), backingfile = NULL){
  G <- x$genotypes
  ind.row <- bigstatsr::rows_along(G)[ind.row]
  ind.col <- bigstatsr::cols_along(G)[ind.col]
  # check_args()
  if (is.null(x$fam)) {
    new_fam <- NULL
  }
  else {
    new_fam <- x$fam[ind.row, , drop = FALSE]
    rownames(new_fam) <- bigstatsr::rows_along(new_fam)
  }
  if (is.null(x$map)) {
    new_map <- NULL
  }
  else {
    new_map <- x$map[ind.col, , drop = FALSE]
    rownames(new_map) <- bigstatsr::rows_along(new_map)
  }
  if (is.null(backingfile)) 
    backingfile <- bigstatsr:::getNewFile(x, "sub")
  G2 <- bigstatsr:::FBM.code256(nrow = length(ind.row), ncol = length(ind.col), 
                                code = G$code256, init = NULL, backingfile = backingfile, 
                                create_bk = TRUE)
  bigsnpr:::replaceSNP(G2, G, rowInd = ind.row, colInd = ind.col)
  snp.list <- structure(list(genotypes = G2, fam = new_fam, 
                             map = new_map), class = "bigSNP")
  rds <- bigstatsr:::sub_bk(G2$backingfile, ".rds")
  saveRDS(snp.list, rds)
  rds
}

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "Asthma"
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_act_Asthma.regenie", sep="")
}else if(trait == 2){
  trait <- "CAD"
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_act_CAD.regenie", sep="")
}else if(trait == 3){
  trait <- "T2D"
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_act_T2D.regenie", sep="")
}else if(trait == 4){
  trait <- "Breast"
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_bp_Breast.regenie", sep="")
}else{
  trait <- "Prostate"
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_bp_Prostate.regenie", sep="")
}

map <- NULL

ldr <- 3/1000
ncores <- 1

colnames(dat) <- c("CHR","POS","SNP_ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
dat$PVAL <- 10^(-1*dat$LOG10P)

write.table(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/",trait,"_assoc.txt"),col.names = T,row.names = F,quote=F)

sumstats <- dat[,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")

obj.bigSNP <- snp_attach("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/all_chr_reference.rds")
NCORES <-  1

for(i in 1:22){
  file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".rds"))
  file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".bk"))
  snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = which(sumstats$chr == i),backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i)) 
}

for(i in 1:22){
  obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".rds"))
  
  map_new <- obj.bigSNP_new$map[-c(3)]
  names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
  
  G   <- obj.bigSNP_new$genotypes
  CHR <- obj.bigSNP_new$map$chromosome
  POS <- obj.bigSNP_new$map$physical.pos
  POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Binary/LDPred2_Genetic_Mappings/", ncores = ncores)
  
  corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
  
  if(anyNA(corr0@x)){
    b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
    b <- as.numeric(names(table(b))[table(b) > 2])
    
    if(length(b) == 0){
      b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
      b <- as.numeric(names(table(b))[table(b) > 1])
    }
    
    file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".rds"))
    file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".bk"))
    
    snp_subset(obj.bigSNP_new,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i))
    
    obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".rds"))
    
    map_new <- obj.bigSNP_new$map[-c(3)]
    names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
    
    G   <- obj.bigSNP_new$genotypes
    CHR <- obj.bigSNP_new$map$chromosome
    POS <- obj.bigSNP_new$map$physical.pos
    POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Binary/LDPred2_Genetic_Mappings/", ncores = ncores)
    
    corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
    while(anyNA(corr0@x)){
      b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
      b <- as.numeric(names(table(b))[table(b) > 1])
      
      file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".rds"))
      file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".bk"))
      
      snp_subset(obj.bigSNP_new,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i))
      
      obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_",trait,"_chr",i,".rds"))
      
      map_new <- obj.bigSNP_new$map[-c(3)]
      names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
      
      G   <- obj.bigSNP_new$genotypes
      CHR <- obj.bigSNP_new$map$chromosome
      POS <- obj.bigSNP_new$map$physical.pos
      POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Binary/LDPred2_Genetic_Mappings/", ncores = ncores)
      
      corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
    }
  }
  
  map <- rbind(map,map_new)
  
  if(i == 1){
    corr <- as_SFBM(corr0, tempfile(), compact = TRUE)
  }else{
    corr$add_columns(corr0, nrow(corr))
  }
}

sumstats <- sumstats[sumstats$rsid %in% map$rsid,]

info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

corr <- corr[info_snp$`_NUM_ID_`,info_snp$`_NUM_ID_`]
corr <- as_SFBM(as(corr, "generalMatrix"))


# Automatic model
ldsc <- snp_ldsc2(corr, df_beta)
h2_est <- ldsc[["h2"]]
print(paste0('Complete data preparation'))

## LDpred2
h2_seq <- seq(0.1,1.5,by = 0.1)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = signif(abs(h2_est) * h2_seq, 3), sparse = c(FALSE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = NCORES)
rownames(beta_grid) = info_snp$rsid
beta_grid = cbind(beta_grid, info_snp[,c('a0','a1','rsid')])
colnames(beta_grid) = c(paste0('e',1:nrow(params)), 'a0','a1','rsid')
beta_grid[is.na(beta_grid)] = 0
beta_grid = as.data.frame(beta_grid)

## LASSOSUM2
delta_path <- function (max=100, min=0.5, n=10){
  sqrt_max <- max^(1/3)
  sqrt_min <- min^(1/3)
  path <- numeric(n)
  for (i in 1:n) {
    path[n+1-i] = (sqrt_max - (sqrt_max-sqrt_min)/(n-1)  * (i-1) )^3;
  }
  return(path)
}
beta_lassosum2 <- snp_lassosum2(corr, df_beta, delta = delta_path(max=100,min=0.5,n=10),ncores = NCORES, maxiter=1000)
params2 <- attr(beta_lassosum2, "grid_param")
rownames(beta_lassosum2) = info_snp$rsid
beta_lassosum2 = cbind(beta_lassosum2, info_snp[,c('a0','a1','rsid')])
colnames(beta_lassosum2) = c(paste0('e',1:nrow(params2)), 'a0','a1','rsid')
beta_lassosum2[is.na(beta_lassosum2)] <- 0
beta_lassosum2 = as.data.frame(beta_lassosum2)

rm(corr0, corr)
print(paste0('Complete'))

# -------- PRS:

## LDpred2 
prs.file <- data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, REF = beta_grid$a1, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2.txt"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_train"))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_tune"))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_validation"))

## LASSOsum2
prs.file <- data.frame(SNP = beta_lassosum2$rsid, ALT = beta_lassosum2$a0, REF = beta_lassosum2$a1, BETA = beta_lassosum2[,1:(ncol(beta_lassosum2)-3)])
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2.txt"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_train"))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_tune"))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_validation"))

################

prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_train.sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_train.log")))
prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_tune.sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_tune.log")))
prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_validation.sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_validation.log")))

sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate AUC for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
AUC_tun_vec <- rep(0,nrow(sets))
if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <- paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}

for(k in 1:nrow(sets)){
  d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",k,"_SUM"))]
  d[,paste0("SCORE",k,"_SUM")] <- -1*d[,paste0("SCORE",k,"_SUM")]
  
  roc_obj <- roc.binary(status = trait,
                       variable = paste0("SCORE",k,"_SUM"),
                       confounders = as.formula(confounders),
                       data = d,
                       precision=seq(0.05,0.95, by=0.05))
  
  AUC_tun_vec[k] <- roc_obj$auc
}

idx <- which.max(AUC_tun_vec)

best_prs_train <- pheno_train[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_tune <- pheno_tuning[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_validation <- pheno_vad[,c("IID",paste0("SCORE",idx,"_SUM"))]

write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)

##### Final Coefficients

all_betas <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2.txt"), sep="")
colnames(all_betas) <- c("SNP","ALT","REF",paste0("LDPred2_SCORE",1:nrow(sets),"_SUM"))
system(paste("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Binary/",trait,"_ldpred2.txt")))

if(trait == "Asthma"){
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_act_Asthma.regenie", sep="")
}else if(trait == "CAD"){
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_act_CAD.regenie", sep="")
}else if(trait == "T2D"){
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_act_T2D.regenie", sep="")
}else if(trait == "Breast"){
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_bp_Breast.regenie", sep="")
}else{
  dat <- read.csv("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_bp_Prostate.regenie", sep="")
}

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")

dat <- dat[,c("CHROM","ID","REF","POS","ALT")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1")


dat <- left_join(dat,all_betas)
dat[is.na(dat)] <- 0

write.csv(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_Final_Coefficients.csv"),row.names = FALSE)





load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_validation_raw <- pheno_vad
pheno_validation_adjusted <- pheno_vad

mod <- lm(as.formula(paste0(paste0("SCORE",idx,"_SUM"),"~pc1 + pc2 + pc3 + pc4 + pc5")),data = pheno_validation_adjusted)
R <- mod$residuals
tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
y_hat <- predict(mod,tmp)
if(sum(sqrt(y_hat)) == 0){
  pheno_validation_adjusted[,paste0("SCORE",idx,"_SUM")] <- 0
}else{
  pheno_validation_adjusted[,paste0("SCORE",idx,"_SUM")] <- R/sqrt(y_hat)
}

pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_raw_NonEUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_raw_UNK <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_MIX <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_NonEUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_adjusted_UNK <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_MIX <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_NonEUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_NonEUR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_UNK[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_UNK[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_MIX[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_MIX[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")])

beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
se_validation_raw_EUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))$coefficients[2,2]
beta_validation_raw_NonEUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_NonEUR,family = binomial()))[2]
se_validation_raw_NonEUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_NonEUR,family = binomial()))$coefficients[2,2]
beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
se_validation_raw_SAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))$coefficients[2,2]
beta_validation_raw_MIX <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_MIX,family = binomial()))[2]
se_validation_raw_MIX <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_MIX,family = binomial()))$coefficients[2,2]
beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
se_validation_raw_AFR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))$coefficients[2,2]
beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
se_validation_raw_EAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))$coefficients[2,2]

beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
se_validation_adjusted_EUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_NonEUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_NonEUR,family = binomial()))[2]
se_validation_adjusted_NonEUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_NonEUR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
se_validation_adjusted_SAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_MIX <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_MIX,family = binomial()))[2]
se_validation_adjusted_MIX <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_MIX,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
se_validation_adjusted_AFR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
se_validation_adjusted_EAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))$coefficients[2,2]

ldpred2_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                              beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_NonEUR,beta_validation_raw_UNK,beta_validation_raw_SAS,beta_validation_raw_MIX,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                              se_raw = c(se_validation_raw_EUR,se_validation_raw_NonEUR,se_validation_raw_UNK,se_validation_raw_SAS,se_validation_raw_MIX,se_validation_raw_AFR,se_validation_raw_EAS), 
                              beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_NonEUR,beta_validation_adjusted_UNK,beta_validation_adjusted_SAS,beta_validation_adjusted_MIX,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                              se_adjusted = c(se_validation_adjusted_EUR,se_validation_adjusted_NonEUR,se_validation_adjusted_UNK,se_validation_adjusted_SAS,se_validation_adjusted_MIX,se_validation_adjusted_AFR,se_validation_adjusted_EAS))

write.csv(ldpred2_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"Best_Betas.csv"),row.names = FALSE)









prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_train.sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_train.log")))
prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_tune.sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_tune.log")))
prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_validation.sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_validation.log")))

## Pull in Phenotypes/Covariates 
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

AUC_tun_vec <- rep(0,300)
for(k in 1:300){
  d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",k,"_SUM"))]
  d[,paste0("SCORE",k,"_SUM")] <- -1*d[,paste0("SCORE",k,"_SUM")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",k,"_SUM"),
                        confounders = as.formula(confounders),
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  
  AUC_tun_vec[k] <- roc_obj$auc
}

idx <- which.max(AUC_tun_vec)

best_prs_train <- pheno_train[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_tune <- pheno_tuning[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_validation <- pheno_vad[,c("IID",paste0("SCORE",idx,"_SUM"))]

write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_lassosum2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad$y_validation <- NA
pheno_vad$y_validation[!is.na(pheno_vad[,trait])] <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad)$residual

pheno_validation_raw <- pheno_vad
pheno_validation_adjusted <- pheno_vad

mod <- lm(as.formula(paste0(paste0("SCORE",idx,"_SUM"),"~pc1 + pc2 + pc3 + pc4 + pc5")),data = pheno_validation_adjusted)
R <- mod$residuals
tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
mod <- lm(y~.,data = tmp)
y_hat <- predict(mod,tmp)
if(sum(sqrt(y_hat)) == 0){
  pheno_validation_adjusted[,paste0("SCORE",idx,"_SUM")] <- 0
}else{
  pheno_validation_adjusted[,paste0("SCORE",idx,"_SUM")] <- R/sqrt(y_hat)
}

pheno_validation_raw_EUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_raw_NonEUR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_raw_UNK <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_MIX <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_NonEUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_validation_adjusted_UNK <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_MIX <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_NonEUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_NonEUR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_UNK[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_UNK[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_MIX[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_MIX[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")])

beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
se_validation_raw_EUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))$coefficients[2,2]
beta_validation_raw_NonEUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_NonEUR,family = binomial()))[2]
se_validation_raw_NonEUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_NonEUR,family = binomial()))$coefficients[2,2]
beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
se_validation_raw_SAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))$coefficients[2,2]
beta_validation_raw_MIX <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_MIX,family = binomial()))[2]
se_validation_raw_MIX <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_MIX,family = binomial()))$coefficients[2,2]
beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
se_validation_raw_AFR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))$coefficients[2,2]
beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
se_validation_raw_EAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))$coefficients[2,2]

beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
se_validation_adjusted_EUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_NonEUR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_NonEUR,family = binomial()))[2]
se_validation_adjusted_NonEUR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_NonEUR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
se_validation_adjusted_SAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_MIX <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_MIX,family = binomial()))[2]
se_validation_adjusted_MIX <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_MIX,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
se_validation_adjusted_AFR <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))$coefficients[2,2]
beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
se_validation_adjusted_EAS <- summary(glm(as.formula(paste0(trait,paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))$coefficients[2,2]

lassosum2_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                                beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_NonEUR,beta_validation_raw_UNK,beta_validation_raw_SAS,beta_validation_raw_MIX,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                se_raw = c(se_validation_raw_EUR,se_validation_raw_NonEUR,se_validation_raw_UNK,se_validation_raw_SAS,se_validation_raw_MIX,se_validation_raw_AFR,se_validation_raw_EAS), 
                                beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_NonEUR,beta_validation_adjusted_UNK,beta_validation_adjusted_SAS,beta_validation_adjusted_MIX,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                se_adjusted = c(se_validation_adjusted_EUR,se_validation_adjusted_NonEUR,se_validation_adjusted_UNK,se_validation_adjusted_SAS,se_validation_adjusted_MIX,se_validation_adjusted_AFR,se_validation_adjusted_EAS))

write.csv(lassosum2_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"Best_Betas.csv"),row.names = FALSE)
