rm(list = ls())

i <- 1

for(i in 1:22){
  rm(list=setdiff(ls(), "i"))
  library(data.table)
  library(dplyr)
  library(pROC)
  library(bigsnpr)
  library(bigsparser)
  library(readr)
  race = c('EUR','AFR','AMR')[2]
  trait = c("HDL","LDL","logTG","TC")[1]
  chr =  i
  ldr = 3/1000
  ncores = 1
  
  temdir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait)
  if (!dir.exists(temdir)){dir.create(temdir)}
  temdir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race)
  if (!dir.exists(temdir)){dir.create(temdir)}
  
  tedir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/intermediate/')
  if (!dir.exists(tedir)){dir.create(tedir)}
  
  temdir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/ldpred2/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  temdir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/lassosum2/')
  if (!dir.exists(temdir)){dir.create(temdir)}
  
  sumraw = as.data.frame(fread(paste0("/data/BB_Bioinformatics/ProjectData/GLGC_cleaned/",race,"/",trait,".txt"),header=T))
  sumraw = sumraw[,c('rsID','CHR','POS_b37','N','BETA','SE','P','A1','A2')]
  colnames(sumraw) = c('SNP_ID','CHR','POS','N','BETA','SE','PVAL','REF','ALT') # REF: allele corresponding to BETA
  
  valdat = bigreadr::fread2(paste0('/data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/',race,'/chr',chr,'.bim'))
  sumraw = sumraw[sumraw$SNP_ID %in% valdat[,2],]
  
  #### read in reference data, this should match as this is what the reference data was in CT
  refdir = paste0('/data/BB_Bioinformatics/ProjectData/GRCh37/',race)
  obj.bigSNP <- snp_attach(paste0(refdir,'/chr',chr,'.rds'))
  map <- obj.bigSNP$map[-c(3)]
  names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
  
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/intermediate/'), ncores = ncores)
  NCORES <-  1
  
  sumstats = sumraw[sumraw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
  set.seed(2020)
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
  sumstats = sumstats[sumstats$rsid %in% map$rsid,]
  info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  rownames(info_snp) = info_snp$rsid
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2,infos.pos = POS2[ind.chr2], size =  ldr)
  corr <- as_SFBM(as(corr0, "dgCMatrix"))
  
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
  
  tedir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/ldpred2/effect/')
  if (!dir.exists(tedir)) dir.create(tedir)
  tedir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/ldpred2/score/')
  if (!dir.exists(tedir)) dir.create(tedir)
  outdir_ldpred2 = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/ldpred2/')
  
  outfile_ldpred2 = paste0(outdir_ldpred2, 'score/ldpred2-chr', chr,'.sscore')
  
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
  
  tedir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/lassosum2/effect/')
  if (!dir.exists(tedir)) dir.create(tedir)
  tedir = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/lassosum2/score/')
  if (!dir.exists(tedir)) dir.create(tedir)
  outdir_lassosum2 = paste0('/data/williamsjacr/multi_ethnic/data/GRCh37/',trait,'/',race,'/lassosum2/')
  
  outfile_lassosum2 = paste0(outdir_lassosum2, 'score/lassosum2-chr', chr,'.sscore')
  
  rm(corr0, corr)
  print(paste0('Complete'))
  
  # -------- PRS:
  
  ## LDpred2 
  prs.file = data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
  tem = bigreadr::fread2(paste0('/data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/',race,'/chr',chr,'.bim'))
  dupid = tem[duplicated(tem[,2]),2]
  prs.file = prs.file[!(prs.file$SNP %in% dupid),]
  write.table(prs.file,file = paste0(outdir_ldpred2,'effect/ldpred2-chr',chr,'.txt'),
              col.names = T,row.names = F,quote=F)
  
  prscode = paste(paste0('/data/williamsjacr/software/plink2'),
                  paste0('--score ',  outdir_ldpred2,'effect/ldpred2-chr',chr,'.txt'),
                  'cols=+scoresums,-scoreavgs header no-mean-imputation ',
                  paste0('--score-col-nums 3-',ncol(prs.file)),
                  paste0(' --bfile /data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/',race,'/chr',chr),
                  " --threads 1",
                  paste0(' --out ', outdir_ldpred2, 'score/ldpred2-chr', chr))
  
  system(prscode)
  
  ## LASSOsum2
  prs.file = data.frame(SNP = beta_lassosum2$rsid, ALT = beta_lassosum2$a0, BETA = beta_lassosum2[,1:(ncol(beta_lassosum2)-3)])
  tem = bigreadr::fread2(paste0('/data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/',race,'/chr',chr,'.bim'))
  dupid = tem[duplicated(tem[,2]),2]
  prs.file = prs.file[!(prs.file$SNP %in% dupid),]
  write.table(prs.file,file = paste0(outdir_lassosum2,'effect/lassosum2-chr',chr,'.txt'),
              col.names = T,row.names = F,quote=F)
  
  prscode = paste(paste0('/data/williamsjacr/software/plink2'),
                  paste0('--score ',  outdir_lassosum2,'effect/lassosum2-chr',chr,'.txt'),
                  'cols=+scoresums,-scoreavgs header no-mean-imputation ',
                  paste0('--score-col-nums 3-',ncol(prs.file)),
                  paste0(' --bfile /data/BB_Bioinformatics/ProjectData/UKBB_sub/genotype/all_data/',race,'/chr',chr),
                  " --threads 1",
                  paste0(' --out ', outdir_lassosum2, 'score/lassosum2-chr', chr))
  
  system(prscode)
  
  print(paste0('Complete Calculating Scores'))
}

#########################################################################################

rm(list = ls())
library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)
for(i in 1:22){
  if(i == 1){
    prs_mat <- fread(paste0("/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/ldpred2/score/ldpred2-chr",i,".sscore"))
  }else{
    prs_mat[,5:55] <- prs_mat[,5:55] + fread(paste0("/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/ldpred2/score/ldpred2-chr",i,".sscore"))[,5:55]
  }
}
# no NA's
# apply(prs_mat,2,function(x){sum(is.na(x))})
rm(i)

write.table(prs_mat,file="/data/williamsjacr/multi_ethnic/data/GRCh37/HDL/AFR/ldpred2/ldpred2_prs_all.txt",sep = "\t")

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
