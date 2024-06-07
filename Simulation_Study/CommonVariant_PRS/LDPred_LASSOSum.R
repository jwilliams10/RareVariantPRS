rm(list=ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)
library(boot)
library(bigstatsr)

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

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

ldr <- 3/1000
ncores <- 1

dat <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/GWAS_Summary_Statistics/Y_Train",i,".Y.glm.linear"), header=FALSE, comment.char="#")
colnames(dat) <- c("CHR","POS","SNP_ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","N","BETA","SE","T_STAT","PVAL","ERRCODE")
dat <- dat[dat$TEST == "ADD",]

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/GWAS_Summary_Statistics/Y_Train",i,".Y.glm.linear")))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/GWAS_Summary_Statistics/Y_Train",i,".log")))

obj.bigSNP <- snp_attach("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.rds")
map <- obj.bigSNP$map[-c(3)]
names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Simulation/Simulation1/LDPred2_Genetic_Mappings/", ncores = ncores)
NCORES <-  1

sumstats <- dat[,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
sumstats <- sumstats[sumstats$rsid %in% map$rsid,]

info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)

if(anyNA(corr0@x)){
  b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
  b <- as.numeric(names(table(b))[table(b) > 2])
  
  if(file.exists(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))){
    file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))
    file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".bk"))
  }
  
  snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i))
  
  obj.bigSNP <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))
  map <- obj.bigSNP$map[-c(3)]
  names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
  
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Simulation/Simulation1/LDPred2_Genetic_Mappings/", ncores = ncores)
  NCORES <-  1
  
  sumstats <- dat[,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
  set.seed(2020)
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
  sumstats <- sumstats[sumstats$rsid %in% map$rsid,]
  
  info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  rownames(info_snp) = info_snp$rsid
  
  df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
  
  corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
  if(anyNA(corr0@x)){
    b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
    b <- as.numeric(names(table(b))[table(b) > 1])
    
    if(file.exists(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))){
      file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))
      file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".bk"))
    }
    
    snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i))
    
    obj.bigSNP <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))
    map <- obj.bigSNP$map[-c(3)]
    names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
    
    G   <- obj.bigSNP$genotypes
    CHR <- obj.bigSNP$map$chromosome
    POS <- obj.bigSNP$map$physical.pos
    POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Simulation/Simulation1/LDPred2_Genetic_Mappings/", ncores = ncores)
    NCORES <-  1
    
    sumstats <- dat[,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
    set.seed(2020)
    names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
    sumstats <- sumstats[sumstats$rsid %in% map$rsid,]
    
    info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
    rownames(info_snp) = info_snp$rsid
    
    df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
    
    corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
  }
}

file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".rds"))
file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i,".bk"))


# Imputing fixes the problem
# G2 = snp_fastImputeSimple(G)
# corr0 <- snp_cor(G2, ind.col = ind.chr2,infos.pos = POS2[ind.chr2], size =  ldr)

corr <- as_SFBM(as(corr0, "generalMatrix"))

# Automatic model
ldsc <- snp_ldsc2(corr0, df_beta)
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
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2",i,".txt"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2",i,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_train",i))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2",i,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_tune",i))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2",i,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_validation",i))

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2",i,".txt")))

## LASSOsum2
prs.file <- data.frame(SNP = beta_lassosum2$rsid, ALT = beta_lassosum2$a0, REF = beta_lassosum2$a1, BETA = beta_lassosum2[,1:(ncol(beta_lassosum2)-3)])
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2",i,".txt"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2",i,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_train",i))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2",i,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_tune",i))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2",i,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common --keep /data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_validation",i))

system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2",i,".txt")))

################

prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_train",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_train",i,".log")))
prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_tune",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_tune",i,".log")))
prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_validation",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/prs_validation",i,".log")))

sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- Y_tune[[i]]
colnames(pheno_tuning) <- c("IID","Y")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,nrow(sets))
model.null <- lm(Y~1,data=pheno_tuning)
for(k in 1:nrow(sets)){
  prs <- pheno_tuning[,paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

best_prs_train <- pheno_train[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_tune <- pheno_tuning[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_validation <- pheno_vad[,c("IID",paste0("SCORE",idx,"_SUM"))]

write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_train_prs_best",i,".txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_tune_prs_best",i,".txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_validation_prs_best",i,".txt"),sep = "\t",row.names = FALSE)






load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad$y_validation <- NA
pheno_vad$y_validation <- lm(Y~1,data=pheno_vad)$residual

pheno_validation_raw <- pheno_vad
pheno_validation_adjusted <- pheno_vad
pheno_validation_adjusted <- inner_join(pheno_validation_adjusted,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])

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
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
pheno_validation_raw_AMR$y_validation <- scale(pheno_validation_raw_AMR$y_validation)
pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)
pheno_validation_raw_EAS$y_validation <- scale(pheno_validation_raw_EAS$y_validation)

pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")])

pheno_validation_adjusted_EUR$y_validation <- scale(pheno_validation_adjusted_EUR$y_validation)
pheno_validation_adjusted_SAS$y_validation <- scale(pheno_validation_adjusted_SAS$y_validation)
pheno_validation_adjusted_AMR$y_validation <- scale(pheno_validation_adjusted_AMR$y_validation)
pheno_validation_adjusted_AFR$y_validation <- scale(pheno_validation_adjusted_AFR$y_validation)
pheno_validation_adjusted_EAS$y_validation <- scale(pheno_validation_adjusted_EAS$y_validation)

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))[2]
se_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))$coefficients[2,2]
beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))[2]
se_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))$coefficients[2,2]
beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))[2]
se_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))$coefficients[2,2]
beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))[2]
se_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))$coefficients[2,2]
beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))[2]
se_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))$coefficients[2,2]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))[2]
se_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))$coefficients[2,2]
beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))[2]
se_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))$coefficients[2,2]
beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))[2]
se_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))$coefficients[2,2]
beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))[2]
se_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))$coefficients[2,2]
beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))[2]
se_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))$coefficients[2,2]

ldpred2_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                              beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                              se_raw = c(se_validation_raw_EUR,se_validation_raw_SAS,se_validation_raw_AMR,se_validation_raw_AFR,se_validation_raw_EAS), 
                              beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                              se_adjusted = c(se_validation_adjusted_EUR,se_validation_adjusted_SAS,se_validation_adjusted_AMR,se_validation_adjusted_AFR,se_validation_adjusted_EAS))

write.csv(ldpred2_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/Best_Betas",i,".csv"),row.names = FALSE)





############################### LASSOSUM2
prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_train",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_train",i,".log")))
prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_tune",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_tune",i,".log")))
prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_validation",i,".sscore"))
system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/prs_validation",i,".log")))

## Pull in Phenotypes/Covariates 
load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- Y_tune[[i]]
colnames(pheno_tuning) <- c("IID","Y")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_vad <- Y_validation[[i]]
colnames(pheno_vad) <- c("IID","Y")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

r2_tun_vec <- rep(0,300)
for(k in 1:300){
  prs <- pheno_tuning[,paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)




load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad$y_validation <- NA
pheno_vad$y_validation <- lm(Y~1,data=pheno_vad)$residual

pheno_validation_raw <- pheno_vad
pheno_validation_adjusted <- pheno_vad
pheno_validation_adjusted <- inner_join(pheno_validation_adjusted,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])

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
pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_validation_raw_EUR$y_validation <- scale(pheno_validation_raw_EUR$y_validation)
pheno_validation_raw_SAS$y_validation <- scale(pheno_validation_raw_SAS$y_validation)
pheno_validation_raw_AMR$y_validation <- scale(pheno_validation_raw_AMR$y_validation)
pheno_validation_raw_AFR$y_validation <- scale(pheno_validation_raw_AFR$y_validation)
pheno_validation_raw_EAS$y_validation <- scale(pheno_validation_raw_EAS$y_validation)

pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")])
pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")])

pheno_validation_adjusted_EUR$y_validation <- scale(pheno_validation_adjusted_EUR$y_validation)
pheno_validation_adjusted_SAS$y_validation <- scale(pheno_validation_adjusted_SAS$y_validation)
pheno_validation_adjusted_AMR$y_validation <- scale(pheno_validation_adjusted_AMR$y_validation)
pheno_validation_adjusted_AFR$y_validation <- scale(pheno_validation_adjusted_AFR$y_validation)
pheno_validation_adjusted_EAS$y_validation <- scale(pheno_validation_adjusted_EAS$y_validation)

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))[2]
se_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))$coefficients[2,2]
beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))[2]
se_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))$coefficients[2,2]
beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))[2]
se_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))$coefficients[2,2]
beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))[2]
se_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))$coefficients[2,2]
beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))[2]
se_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))$coefficients[2,2]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))[2]
se_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))$coefficients[2,2]
beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))[2]
se_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))$coefficients[2,2]
beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))[2]
se_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))$coefficients[2,2]
beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))[2]
se_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))$coefficients[2,2]
beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))[2]
se_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))$coefficients[2,2]

lassosum2_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                se_raw = c(se_validation_raw_EUR,se_validation_raw_SAS,se_validation_raw_AMR,se_validation_raw_AFR,se_validation_raw_EAS), 
                                beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                se_adjusted = c(se_validation_adjusted_EUR,se_validation_adjusted_SAS,se_validation_adjusted_AMR,se_validation_adjusted_AFR,se_validation_adjusted_EAS))

write.csv(lassosum2_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/Best_Betas",i,".csv"),row.names = FALSE)
