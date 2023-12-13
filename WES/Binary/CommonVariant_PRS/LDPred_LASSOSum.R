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
  if(!file.exists(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_chr",i,".rds")) & !file.exists(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_chr",i,".bk"))){
    file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_chr",i,".rds"))
    file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_chr",i,".bk"))
    snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = which(sumstats$chr == i),backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_chr",i)) 
  }
  
  obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/BEDFiles/reference_chr",i,".rds"))
  
  map_new <- obj.bigSNP_new$map[-c(3)]
  names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
  map <- rbind(map,map_new)
  
  
  G   <- obj.bigSNP_new$genotypes
  CHR <- obj.bigSNP_new$map$chromosome
  POS <- obj.bigSNP_new$map$physical.pos
  POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Binary/LDPred2_Genetic_Mappings/", ncores = ncores)
  
  if(i == 1){
    corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
    corr <- as_SFBM(corr0, tempfile(), compact = TRUE)
  }else{
    corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
    corr$add_columns(corr0, nrow(corr))
  }
}

sumstats <- sumstats[sumstats$rsid %in% map$rsid,]

info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

corr <- corr[info_snp$`_NUM_ID_`,info_snp$`_NUM_ID_`]
corr <- as_SFBM(as(corr, "generalMatrix"))



# if(anyNA(corr0@x)){
#   b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)  
#   # b_list <- NULL
#   # for(j in 1:nrow(b)){
#   #   b_list <- c(b_list,b[j,])
#   # }
#   # b <- unique(b_list)
#   # rm(b_list)
#   b <- as.numeric(names(table(b))[table(b) > 2])
#   
#   b <- info_snp$`_NUM_ID_`[b]
#   
#   if(file.exists(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.rds"))){
#     file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.rds"))
#     file.remove(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.bk"))
#   }
#   
#   snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference",i))
#   
#   obj.bigSNP <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/reference.rds"))
#   map <- obj.bigSNP$map[-c(3)]
#   names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
#   
#   G   <- obj.bigSNP$genotypes
#   CHR <- obj.bigSNP$map$chromosome
#   POS <- obj.bigSNP$map$physical.pos
#   POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Simulation/Simulation1/LDPred2_Genetic_Mappings/", ncores = ncores)
#   NCORES <-  1
#   
#   sumstats <- dat[,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
#   set.seed(2020)
#   names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
#   sumstats <- sumstats[sumstats$rsid %in% map$rsid,]
#   
#   info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
#   rownames(info_snp) = info_snp$rsid
#   
#   df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
#   
#   corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
# }


# Imputing fixes the problem
# G2 = snp_fastImputeSimple(G)
# corr0 <- snp_cor(G2, ind.col = ind.chAUC,infos.pos = POS2[ind.chAUC], size =  ldr)

# corr <- as_SFBM(as(corr, "generalMatrix"))

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
prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_tune.sscore"))
prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_prs_validation.sscore"))

sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

#calculate AUC for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
AUC_tun_vec <- rep(0,nrow(sets))
if(trait %in% c("Breast","Prostate")){
  confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
}else{
  confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
}

for(k in 1:nrow(sets)){
  d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",k,"_SUM"))]
  d[,paste0("SCORE",k,"_SUM")] <- -1*d[,paste0("SCORE",k,"_SUM")]
  
  roc_obj <- roc.binary(status = trait,
                       variable = paste0("SCORE",k,"_SUM"),
                       confounders = confounders,
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

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                     variable = paste0("SCORE",idx,"_SUM"),
                     confounders = confounders,
                     data = d,
                     precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                       variable = paste0("SCORE",idx,"_SUM"),
                       confounders = confounders,
                       data = d_sub,
                       precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LDPred2_EUR",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ldpred2.result <- list(AUC.result,AUC_tun_vec)

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_EUR.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LDPred2_NonEur",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ldpred2.result <- list(AUC.result,AUC_tun_vec)

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_NonEur.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LDPred2_UNK",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ldpred2.result <- list(AUC.result,AUC_tun_vec)

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_UNK.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LDPred2_AFR",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ldpred2.result <- list(AUC.result,AUC_tun_vec)

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_AFR.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LDPred2_SAS",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ldpred2.result <- list(AUC.result,AUC_tun_vec)

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_SAS.RData"))

## bootstrap the AUC to provide an approximate distribution 
if(trait %in% c("Prostate","CAD")){
  ldpred2.result <- NA
}else{
  d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
  d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SCORE",idx,"_SUM"),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  
  AUC.result <- data.frame(method = "LDPred2_EAS",
                           AUC = AUC,
                           AUC_low = ci_result$percent[4],
                           AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ldpred2.result <- list(AUC.result,AUC_tun_vec)
}

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_EAS.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LDPred2_MIX",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
ldpred2.result <- list(AUC.result,AUC_tun_vec)

save(ldpred2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",trait,"_ldpred2_result_MIX.RData"))










prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_train.sscore"))
prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_tune.sscore"))
prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_prs_validation.sscore"))

## Pull in Phenotypes/Covariates 
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

AUC_tun_vec <- rep(0,300)
for(k in 1:300){
  d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",k,"_SUM"))]
  d[,paste0("SCORE",k,"_SUM")] <- -1*d[,paste0("SCORE",k,"_SUM")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",k,"_SUM"),
                        confounders = confounders,
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

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LASSOSUM2_EUR",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EUR.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LASSOSUM2_NonEur",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_NonEur.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LASSOSUM2_UNK",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_UNK.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LASSOSUM2_AFR",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_AFR.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LASSOSUM2_SAS",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_SAS.RData"))

## bootstrap the AUC to provide an approximate distribution 
if(trait %in% c("Prostate","CAD")){
  LASSOSUM2.result <- NA
}else{
  d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
  d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = paste0("SCORE",idx,"_SUM"),
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  
  AUC.result <- data.frame(method = "LASSOSUM2_EAS",
                           AUC = AUC,
                           AUC_low = ci_result$percent[4],
                           AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)
}

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EAS.RData"))

## bootstrap the AUC to provide an approximate distribution 
d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",idx,"_SUM"))]
d[,paste0("SCORE",idx,"_SUM")] <- -1*d[,paste0("SCORE",idx,"_SUM")]

roc_obj <- roc.binary(status = trait,
                      variable = paste0("SCORE",idx,"_SUM"),
                      confounders = confounders,
                      data = d,
                      precision=seq(0.05,0.95, by=0.05))
AUC <- roc_obj$auc

calc_auc <- function(data, indices) {
  d_sub <- data[indices,] # allows boot to select sample
  roc_obj <- roc.binary(status = trait,
                        variable = paste0("SCORE",idx,"_SUM"),
                        confounders = confounders,
                        data = d_sub,
                        precision=seq(0.05,0.95, by=0.05))
  return(roc_obj$auc)
}
boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
ci_result <- boot.ci(boot_AUC, type = "perc")

AUC.result <- data.frame(method = "LASSOSUM2_MIX",
                        AUC = AUC,
                        AUC_low = ci_result$percent[4],
                        AUC_high = ci_result$percent[5]
)

## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
LASSOSUM2.result <- list(AUC.result,AUC_tun_vec)

save(LASSOSUM2.result,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_MIX.RData"))

