rm(list=ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)

chr <- as.numeric(commandArgs(TRUE)[1])
# chr <- 19
ldr <- 3/1000
ncores <- 1

dat <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/all_chr_assoc.LDLadj.norm.glm.linear", header=FALSE, comment.char="#")
colnames(dat) <- c("CHR","POS","SNP_ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","N","BETA","SE","T_STAT","PVAL","ERRCODE")
dat <- dat[dat$TEST == "ADD",]

if(file.exists(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".rds"))){
  file.remove(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".bk"))
  file.remove(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".rds"))
}

#### read in reference data, this should match as this is what the reference data was in CT
if(!file.exists(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".rds"))){
  snp_readBed(paste0("/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",chr,"_common_reference.bed"),backingfile = paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr))
}
obj.bigSNP <- snp_attach(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".rds"))
map <- obj.bigSNP$map[-c(3)]
names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0("/data/williamsjacr/UKB_WES_lipids/Data/LDPred2_Genetic_Mappings/chr",chr), ncores = ncores)
NCORES <-  1

sumstats <- dat[dat$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
set.seed(2020)
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
sumstats <- sumstats[sumstats$rsid %in% map$rsid,]

info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid

ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

corr0 <- snp_cor(G, ind.col = ind.chr2,infos.pos = POS2[ind.chr2], size =  ldr)

# b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)  
# b_list <- NULL
# for(i in 1:nrow(b)){
#   b_list <- c(b_list,b[,i])
# }
# b <- unique(b_list)
# rm(b_list)
# 
# a <- G[1:3000,ind.chr2[b]]
# sum(a[,1] *a[,2],na.rm = TRUE)

if(anyNA(corr0@x)){
  b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)  
  b_list <- NULL
  for(i in 1:nrow(b)){
    b_list <- c(b_list,b[,i])
  }
  b <- unique(b_list)
  rm(b_list)
  
  b <- ind.chr2[b]
    
  file.remove(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".bk"))
  file.remove(paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr,".rds"))
  bigsnpr::snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_lipids/Data/rdsfiles_LDPred/reference_chr",chr))
  
  sumstats <- dat[dat$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
  set.seed(2020)
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
  sumstats <- sumstats[sumstats$rsid %in% map$rsid,]
  
  sumstats <- sumstats[!(sumstats$pos %in% POS[b]),]
  
  info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  rownames(info_snp) = info_snp$rsid
  
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2,infos.pos = POS2[ind.chr2], size =  ldr)
}


# Imputing fixes the problem
# G2 = snp_fastImputeSimple(G)
# corr0 <- snp_cor(G2, ind.col = ind.chr2,infos.pos = POS2[ind.chr2], size =  ldr)

corr <- as_SFBM(as(corr0, "generalMatrix"))

# Automatic model
ldsc <- snp_ldsc2(corr0, df_beta)
h2_est <- ldsc[["h2"]]
print(paste0('Complete data preparation'))

## LDpred2
h2_seq <- c(0.7, 1, 1.4)
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
prs.file <- data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2-chr",chr,".txt"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2-chr",chr,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_train --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_train_chr",chr))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2-chr",chr,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_tune --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_tune_chr",chr))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2-chr",chr,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_validation --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/prs_validation_chr",chr))

## LASSOsum2
prs.file <- data.frame(SNP = beta_lassosum2$rsid, ALT = beta_lassosum2$a0, BETA = beta_lassosum2[,1:(ncol(beta_lassosum2)-3)])
write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2-chr",chr,".txt"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2-chr",chr,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_train --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_train_chr",chr))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2-chr",chr,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_tune --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_tune_chr",chr))
system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/lassosum2-chr",chr,".txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 3-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_validation --threads 1 --out /data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/prs_validation_chr",chr))
