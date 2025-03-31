rm(list=ls())

# for array in 1 2 3 4 5 6;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/CommonVariant_PRS/LDPred_LASSOSum.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/CommonVariant_PRS/LDPred_LASSOSum.sh -icmd="bash LDPred_LASSOSum.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/ --priority high --instance-type mem3_ssd1_v2_x16 
# done

if(!("data.table" %in% rownames(installed.packages()))){
  install.packages("data.table",quiet = TRUE)
}

if(!("dplyr" %in% rownames(installed.packages()))){
  install.packages("dplyr",quiet = TRUE)
}

if(!("bigsparser" %in% rownames(installed.packages()))){
  install.packages("bigsparser",quiet = TRUE)
}

if(!("readr" %in% rownames(installed.packages()))){
  install.packages("readr",quiet = TRUE)
}

if(!("boot" %in% rownames(installed.packages()))){
  install.packages("boot",quiet = TRUE)
}

if(!("RISCA" %in% rownames(installed.packages()))){
  install.packages("RISCA",quiet = TRUE)
}

if(!("pROC" %in% rownames(installed.packages()))){
  install.packages("pROC",quiet = TRUE)
}

if(!("remotes" %in% rownames(installed.packages()))){
  install.packages("remotes",quiet = TRUE)
}

if(!("bigstatsr" %in% rownames(installed.packages()))){
  install.packages("bigstatsr",quiet = TRUE)
}

if(!("bigsnpr" %in% rownames(installed.packages()))){
  remotes::install_github("privefl/bigsnpr")
}

if(!("R.utils" %in% rownames(installed.packages()))){
  install.packages("R.utils",quiet = TRUE)
}

if(!("stringr" %in% rownames(installed.packages()))){
  install.packages("stringr",quiet = TRUE)
}

library(data.table)
library(dplyr)
library(pROC)
library(bigsnpr)
library(bigsparser)
library(readr)
library(boot)
library(bigstatsr)
library(stringr)

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
  trait <- "BMI"
  dat <- read.csv("regenie_step2_continuous_BMI.regenie", sep="")
}else if(trait == 2){
  trait <- "TC"
  dat <- read.csv("regenie_step2_continuous_TC.regenie", sep="")
}else if(trait == 3){
  trait <- "HDL"
  dat <- read.csv("regenie_step2_continuous_HDL.regenie", sep="")
}else if(trait == 4){
  trait <- "LDL"
  dat <- read.csv("regenie_step2_continuous_LDL.regenie", sep="")
}else if(trait == 5){
  trait <- "logTG"
  dat <- read.csv("regenie_step2_continuous_logTG.regenie", sep="")
}else{
  trait <- "Height"
  dat <- read.csv("regenie_step2_continuous_Height.regenie", sep="")
}

map <- NULL

ldr <- 3/1000
ncores <- 1

colnames(dat) <- c("CHR","POS","SNP_ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
dat$PVAL <- 10^(-1*dat$LOG10P)

dat <- as.data.frame(dat)

sumstats <- dat[,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N',"SNP_ID")]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff","SNP_ID")
rm(dat)

set.seed(2020)

w_hm3 <- read.delim("w_hm3.snplist")

system("rm w_hm3.snplist")

SNP_GRCh37_38_match_update <- readRDS("SNP_GRCh37_38_match_update.rds")

system("rm SNP_GRCh37_38_match_update.rds")

unique_id <- paste0(SNP_GRCh37_38_match_update$chr,"_",SNP_GRCh37_38_match_update$pos38,"_",
                    toupper(SNP_GRCh37_38_match_update$allele1_38),"_",toupper(SNP_GRCh37_38_match_update$allele2_38))

merge_data <- data.frame(rsid = SNP_GRCh37_38_match_update$rsid, unique_id = unique_id)

rm(SNP_GRCh37_38_match_update);rm(unique_id)

sumstats$CHR_POS <- paste0(sumstats$chr,"_",sumstats$pos)
sumstats$unique_id1 <- paste0(sumstats$CHR_POS,"_",
                              toupper(sumstats$a0),"_",toupper(sumstats$a1))
sumstats$unique_id2 <- paste0(sumstats$CHR_POS,"_",
                              toupper(sumstats$a1),"_",toupper(sumstats$a0))

sumstats <- subset(sumstats,select = -rsid)
sumstats$SNPID <- NA

sumstats <- left_join(sumstats,merge_data,by = c("unique_id1" = "unique_id"))
sumstats$SNPID[!is.na(sumstats$rsid)] <- sumstats$rsid[!is.na(sumstats$rsid)]

sumstats <- subset(sumstats,select = -c(rsid))

sumstats <- left_join(sumstats,merge_data,by = c("unique_id2" = "unique_id"))
sumstats$SNPID[!is.na(sumstats$rsid)] <- sumstats$rsid[!is.na(sumstats$rsid)]

rm(merge_data)

gc()

sumstats <- sumstats[,c('chr', 'SNPID', 'pos', 'a0', 'a1', 'beta', 'beta_se', 'p', 'n_eff',"SNP_ID")]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff","SNP_ID")

idx <- sumstats$rsid %in% w_hm3$SNP

subsetted_sumstats <- sumstats[sumstats$rsid %in% w_hm3$SNP,]

# system("dx download UKB_PRS:JW/UKB_Phenotypes/Data/all_chr_reference.rds")
# system("dx download UKB_PRS:JW/UKB_Phenotypes/Data/all_chr_reference.bk")

# system("dx download UKB_PRS:JW/UKB_Phenotypes/Data/LDPred2_Genetic_Mappings.zip")
system("unzip LDPred2_Genetic_Mappings.zip")
system("rm LDPred2_Genetic_Mappings.zip")

obj.bigSNP <- snp_attach("all_chr_reference.rds")
NCORES <-  1

for(i in 1:22){
  snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = which((sumstats$chr == i) & idx),backingfile = paste0("reference_chr",i)) 
}

file.remove("all_chr_reference.rds")
file.remove("all_chr_reference.bk")

for(i in 1:22){
  print(i)
  obj.bigSNP_new <- snp_attach(paste0("reference_chr",i,".rds"))
  
  map_new <- obj.bigSNP_new$map[-c(3)]
  names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
  
  print("extract genotypes")
  G   <- obj.bigSNP_new$genotypes
  CHR <- obj.bigSNP_new$map$chromosome
  POS <- obj.bigSNP_new$map$physical.pos
  print("get positions")
  POS2 <- snp_asGeneticPos(CHR, POS, dir ="LDPred2_Genetic_Mappings", ncores = ncores)
  
  print("build cor")
  corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
  
  if(anyNA(corr0@x)){
    b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
    b <- as.numeric(names(table(b))[table(b) > 2])
    
    if(length(b) == 0){
      b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
      b <- as.numeric(names(table(b))[table(b) > 1])
    }
    
    file.remove(paste0("reference_chr",i,".rds"))
    file.remove(paste0("reference_chr",i,".bk"))
    
    snp_subset(obj.bigSNP_new,ind.row = 1:3000,ind.col = -b,backingfile = paste0("reference_chr",i))
    
    obj.bigSNP_new <- snp_attach(paste0("reference_chr",i,".rds"))
    
    map_new <- obj.bigSNP_new$map[-c(3)]
    names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
    
    G   <- obj.bigSNP_new$genotypes
    CHR <- obj.bigSNP_new$map$chromosome
    POS <- obj.bigSNP_new$map$physical.pos
    POS2 <- snp_asGeneticPos(CHR, POS, dir ="LDPred2_Genetic_Mappings", ncores = ncores)
    
    corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
    while(anyNA(corr0@x)){
      b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
      b <- as.numeric(names(table(b))[table(b) > 1])
      
      file.remove(paste0("reference_chr",i,".rds"))
      file.remove(paste0("reference_chr",i,".bk"))
      
      snp_subset(obj.bigSNP_new,ind.row = 1:3000,ind.col = -b,backingfile = paste0("reference_chr",i))
      
      obj.bigSNP_new <- snp_attach(paste0("reference_chr",i,".rds"))
      
      map_new <- obj.bigSNP_new$map[-c(3)]
      names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
      
      G   <- obj.bigSNP_new$genotypes
      CHR <- obj.bigSNP_new$map$chromosome
      POS <- obj.bigSNP_new$map$physical.pos
      POS2 <- snp_asGeneticPos(CHR, POS, dir ="LDPred2_Genetic_Mappings", ncores = ncores)
      
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

for(i in 1:22){
  file.remove(paste0("reference_chr",i,".rds"))
  file.remove(paste0("reference_chr",i,".bk"))
}

subsetted_sumstats <- subsetted_sumstats[subsetted_sumstats$SNP_ID %in% map$rsid,]

map$SNP_ID <- map$rsid
map <- subset(map,select = -rsid)
map <- inner_join(map,subsetted_sumstats[,c("SNP_ID","rsid")])
map <- subset(map, select = -SNP_ID)

info_snp <- snp_match(subsetted_sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
rownames(info_snp) = info_snp$rsid

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

corr <- corr[info_snp$`_NUM_ID_`,info_snp$`_NUM_ID_`]
corr <- as_SFBM(as(corr, "generalMatrix"))


gc()


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

system("rm -r LDPred2_Genetic_Mappings/")
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


# -------- PRS:

## LDpred2

prs.file <- data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, REF = beta_grid$a1, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
prs.file <- inner_join(prs.file,subsetted_sumstats[,c("SNP_ID","rsid")],by = c("SNP" = "rsid"))
prs.file <- data.frame(SNP = prs.file$SNP_ID, ALT = prs.file$ALT, REF = prs.file$REF, BETA = prs.file[,str_detect(colnames(prs.file),"BETA")])
write.table(prs.file,file = paste0(trait,"_ldpred2.txt"),col.names = T,row.names = F,quote=F)

system(paste0("plink2 --score ",trait,"_ldpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile all_chr --keep train.txt --threads 1 --out ",trait,"_prs_train_ldpred2"))
system(paste0("plink2 --score ",trait,"_ldpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile all_chr --keep tune.txt --threads 1 --out ",trait,"_prs_tune_ldpred2"))
system(paste0("plink2 --score ",trait,"_ldpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile all_chr --keep validation.txt --threads 1 --out ",trait,"_prs_validation_ldpred2"))

## LASSOsum2

prs.file <- data.frame(SNP = beta_lassosum2$rsid, ALT = beta_lassosum2$a0, REF = beta_lassosum2$a1, BETA = beta_lassosum2[,1:(ncol(beta_lassosum2)-3)])
prs.file <- inner_join(prs.file,subsetted_sumstats[,c("SNP_ID","rsid")],by = c("SNP" = "rsid"))
prs.file <- data.frame(SNP = prs.file$SNP_ID, ALT = prs.file$ALT, REF = prs.file$REF, BETA = prs.file[,str_detect(colnames(prs.file),"BETA")])
write.table(prs.file,file = paste0(trait,"_lassosum2.txt"),col.names = T,row.names = F,quote=F)

system(paste0("plink2 --score ",trait,"_lassosum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile all_chr --keep train.txt --threads 1 --out ",trait,"_prs_train_lassosum2"))
system(paste0("plink2 --score ",trait,"_lassosum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile all_chr --keep tune.txt --threads 1 --out ",trait,"_prs_tune_lassosum2"))
system(paste0("plink2 --score ",trait,"_lassosum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile all_chr --keep validation.txt --threads 1 --out ",trait,"_prs_validation_lassosum2"))

################

system(paste0("rm all_chr.bed"))
system(paste0("rm all_chr.bim"))
system(paste0("rm all_chr.fam"))

system(paste0("rm train.txt"))
system(paste0("rm tune.txt"))
system(paste0("rm validation.txt"))


prs_mat_train <- read.delim(paste0("",trait,"_prs_train_ldpred2.sscore"))
prs_mat_tune <- read.delim(paste0("",trait,"_prs_tune_ldpred2.sscore"))
prs_mat_validation <- read.delim(paste0("",trait,"_prs_validation_ldpred2.sscore"))

sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))


pheno_train <- read.delim("All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

#calculate R2 for each of the tuning dataset
# This is done by regressing the residuals of the model with all covariates against the prs
r2_tun_vec <- rep(0,nrow(sets))
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning)
for(k in 1:nrow(sets)){
  prs <- pheno_tuning[!is.na(pheno_tuning[,trait]),paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

best_prs_train <- pheno_train[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_tune <- pheno_tuning[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_validation <- pheno_vad[,c("IID",paste0("SCORE",idx,"_SUM"))]

write.table(best_prs_train,file=paste0("",trait,"_ldpred2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_tune,file=paste0("",trait,"_ldpred2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_validation,file=paste0("",trait,"_ldpred2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)

##### Final Coefficients

all_betas <- read.csv(paste0(trait,"_ldpred2.txt"), sep="")
colnames(all_betas) <- c("SNP","ALT","REF",paste0("LDPred2_SCORE",1:nrow(sets),"_SUM"))
system(paste("rm ",paste0(trait,"_ldpred2.txt")))

if(trait == "BMI"){
  dat <- read.csv("regenie_step2_continuous_BMI.regenie", sep="")
}else if(trait == "TC"){
  dat <- read.csv("regenie_step2_continuous_TC.regenie", sep="")
}else if(trait == "HDL"){
  dat <- read.csv("regenie_step2_continuous_HDL.regenie", sep="")
}else if(trait == "LDL"){
  dat <- read.csv("regenie_step2_continuous_LDL.regenie", sep="")
}else if(trait == "logTG"){
  dat <- read.csv("regenie_step2_continuous_logTG.regenie", sep="")
}else{
  dat <- read.csv("regenie_step2_continuous_Height.regenie", sep="")
}

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")

dat <- dat[,c("CHROM","ID","REF","POS","ALT")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1")

dat <- left_join(dat,all_betas)
dat[is.na(dat)] <- 0

rm(all_betas)
gc()

write.csv(dat,file = paste0(trait,"_Final_Coefficients_LDPred2.csv"),row.names = FALSE)







load("all_phenotypes.RData")

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

Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = boot_data))[2]
  return(c(result))
}

R2_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = boot_data))$r.squared
  return(c(result))
}

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EUR <- sd(boot_R2$t)
R2_lower_validation_raw_EUR <- R2_ci$basic[4]
R2_upper_validation_raw_EUR <- R2_ci$basic[5]

beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_SAS <- sd(boot_R2$t)
R2_lower_validation_raw_SAS <- R2_ci$basic[4]
R2_upper_validation_raw_SAS <- R2_ci$basic[5]

beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AMR <- sd(boot_R2$t)
R2_lower_validation_raw_AMR <- R2_ci$basic[4]
R2_upper_validation_raw_AMR <- R2_ci$basic[5]

beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AFR <- sd(boot_R2$t)
R2_lower_validation_raw_AFR <- R2_ci$basic[4]
R2_upper_validation_raw_AFR <- R2_ci$basic[5]

beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))[2]
boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EAS <- sd(boot_beta$t)
beta_lower_validation_raw_EAS <- beta_ci$basic[4]
beta_upper_validation_raw_EAS <- beta_ci$basic[5]

R2_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EAS <- sd(boot_R2$t)
R2_lower_validation_raw_EAS <- R2_ci$basic[4]
R2_upper_validation_raw_EAS <- R2_ci$basic[5]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
R2_lower_validation_adjusted_EUR <- R2_ci$basic[4]
R2_upper_validation_adjusted_EUR <- R2_ci$basic[5]

beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_SAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_SAS <- R2_ci$basic[5]

beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AMR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AMR <- R2_ci$basic[5]

beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AFR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AFR <- R2_ci$basic[5]

beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]

R2_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_EAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_EAS <- R2_ci$basic[5]

ldpred2_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                              beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                              beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                              beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR,beta_lower_validation_raw_EAS), 
                              beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR,beta_upper_validation_raw_EAS), 
                              R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                              R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                              R2_lower_raw = c(R2_lower_validation_raw_EUR,R2_lower_validation_raw_SAS,R2_lower_validation_raw_AMR,R2_lower_validation_raw_AFR,R2_lower_validation_raw_EAS),
                              R2_upper_raw = c(R2_upper_validation_raw_EUR,R2_upper_validation_raw_SAS,R2_upper_validation_raw_AMR,R2_upper_validation_raw_AFR,R2_upper_validation_raw_EAS),
                              beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                              beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                              beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR,beta_lower_validation_adjusted_EAS), 
                              beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR,beta_upper_validation_adjusted_EAS), 
                              R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                              R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS),
                              R2_lower_adjusted = c(R2_lower_validation_adjusted_EUR,R2_lower_validation_adjusted_SAS,R2_lower_validation_adjusted_AMR,R2_lower_validation_adjusted_AFR,R2_lower_validation_adjusted_EAS),
                              R2_upper_adjusted = c(R2_upper_validation_adjusted_EUR,R2_upper_validation_adjusted_SAS,R2_upper_validation_adjusted_AMR,R2_upper_validation_adjusted_AFR,R2_upper_validation_adjusted_EAS))

write.csv(ldpred2_Results,file = paste0(trait,"Best_Betas_LDPred2.csv"),row.names = FALSE)












prs_mat_train <- read.delim(paste0("",trait,"_prs_train_lassosum2.sscore"))
prs_mat_tune <- read.delim(paste0("",trait,"_prs_tune_lassosum2.sscore"))
prs_mat_validation <- read.delim(paste0("",trait,"_prs_validation_lassosum2.sscore"))

## Pull in Phenotypes/Covariates
pheno_train <- read.delim("All_Train.txt")
pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")

pheno_tuning <- read.delim("All_Tune.txt")
pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")

pheno_vad <- read.delim("All_Validation.txt")
pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")

r2_tun_vec <- rep(0,300)
for(k in 1:300){
  prs <- pheno_tuning[!is.na(pheno_tuning[,trait]),paste0("SCORE",k,"_SUM")]
  model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
  r2_tun_vec[k] <- summary(model.prs)$r.square
}

idx <- which.max(r2_tun_vec)

best_prs_train <- pheno_train[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_tune <- pheno_tuning[,c("IID",paste0("SCORE",idx,"_SUM"))]
best_prs_validation <- pheno_vad[,c("IID",paste0("SCORE",idx,"_SUM"))]

write.table(best_prs_train,file=paste0("",trait,"_lassosum2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_tune,file=paste0("",trait,"_lassosum2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
write.table(best_prs_validation,file=paste0("",trait,"_lassosum2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)

##### Final Coefficients

all_betas <- read.csv(paste0(trait,"_lassosum2.txt"), sep="")
colnames(all_betas) <- c("SNP","ALT","REF",paste0("LASSOSum2_SCORE",1:300,"_SUM"))
system(paste("rm ",paste0(trait,"_lassosum2.txt")))

if(trait == "BMI"){
  dat <- read.csv("regenie_step2_continuous_BMI.regenie", sep="")
  system("rm regenie_step2_continuous_BMI.regenie")
}else if(trait == "TC"){
  dat <- read.csv("regenie_step2_continuous_TC.regenie", sep="")
  system("rm regenie_step2_continuous_TC.regenie")
}else if(trait == "HDL"){
  dat <- read.csv("regenie_step2_continuous_HDL.regenie", sep="")
  system("rm regenie_step2_continuous_HDL.regenie")
}else if(trait == "LDL"){
  dat <- read.csv("regenie_step2_continuous_LDL.regenie", sep="")
  system("rm regenie_step2_continuous_LDL.regenie")
}else if(trait == "logTG"){
  dat <- read.csv("regenie_step2_continuous_logTG.regenie", sep="")
  system("rm regenie_step2_continuous_logTG.regenie")
}else{
  dat <- read.csv("regenie_step2_continuous_Height.regenie", sep="")
  system("rm regenie_step2_continuous_Height.regenie")
}

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")

dat <- dat[,c("CHROM","ID","REF","POS","ALT")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1")

dat <- left_join(dat,all_betas)
dat[is.na(dat)] <- 0

rm(all_betas)
gc()

write.csv(dat,file = paste0(trait,"_Final_Coefficients_LASSOSum.csv"),row.names = FALSE)



load("all_phenotypes.RData")
file.remove("all_phenotypes.RData")

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

Beta_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = boot_data))[2]
  return(c(result))
}

R2_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = boot_data))$r.squared
  return(c(result))
}

beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))[2]
boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EUR <- sd(boot_R2$t)
R2_lower_validation_raw_EUR <- R2_ci$basic[4]
R2_upper_validation_raw_EUR <- R2_ci$basic[5]

beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))[2]
boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_SAS <- sd(boot_R2$t)
R2_lower_validation_raw_SAS <- R2_ci$basic[4]
R2_upper_validation_raw_SAS <- R2_ci$basic[5]

beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))[2]
boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AMR <- sd(boot_R2$t)
R2_lower_validation_raw_AMR <- R2_ci$basic[4]
R2_upper_validation_raw_AMR <- R2_ci$basic[5]

beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))[2]
boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AFR <- sd(boot_R2$t)
R2_lower_validation_raw_AFR <- R2_ci$basic[4]
R2_upper_validation_raw_AFR <- R2_ci$basic[5]

beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))[2]
boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EAS <- sd(boot_beta$t)
beta_lower_validation_raw_EAS <- beta_ci$basic[4]
beta_upper_validation_raw_EAS <- beta_ci$basic[5]

R2_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_raw_EAS))$r.squared
boot_R2 <- boot(data = pheno_validation_raw_EAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EAS <- sd(boot_R2$t)
R2_lower_validation_raw_EAS <- R2_ci$basic[4]
R2_upper_validation_raw_EAS <- R2_ci$basic[5]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EUR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
R2_lower_validation_adjusted_EUR <- R2_ci$basic[4]
R2_upper_validation_adjusted_EUR <- R2_ci$basic[5]

beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_SAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_SAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_SAS <- R2_ci$basic[5]

beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AMR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AMR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AMR <- R2_ci$basic[5]

beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))[2]
boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_AFR))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AFR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AFR <- R2_ci$basic[5]

beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))[2]
boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]

R2_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("SCORE",idx,"_SUM"))),data = pheno_validation_adjusted_EAS))$r.squared
boot_R2 <- boot(data = pheno_validation_adjusted_EAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_EAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_EAS <- R2_ci$basic[5]

lassosum2_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                                beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR,beta_lower_validation_raw_EAS), 
                                beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR,beta_upper_validation_raw_EAS), 
                                R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                                R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                                R2_lower_raw = c(R2_lower_validation_raw_EUR,R2_lower_validation_raw_SAS,R2_lower_validation_raw_AMR,R2_lower_validation_raw_AFR,R2_lower_validation_raw_EAS),
                                R2_upper_raw = c(R2_upper_validation_raw_EUR,R2_upper_validation_raw_SAS,R2_upper_validation_raw_AMR,R2_upper_validation_raw_AFR,R2_upper_validation_raw_EAS),
                                beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                                beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR,beta_lower_validation_adjusted_EAS), 
                                beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR,beta_upper_validation_adjusted_EAS), 
                                R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                                R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS),
                                R2_lower_adjusted = c(R2_lower_validation_adjusted_EUR,R2_lower_validation_adjusted_SAS,R2_lower_validation_adjusted_AMR,R2_lower_validation_adjusted_AFR,R2_lower_validation_adjusted_EAS),
                                R2_upper_adjusted = c(R2_upper_validation_adjusted_EUR,R2_upper_validation_adjusted_SAS,R2_upper_validation_adjusted_AMR,R2_upper_validation_adjusted_AFR,R2_upper_validation_adjusted_EAS))

write.csv(lassosum2_Results,file = paste0(trait,"Best_Betas_LASSOSum.csv"),row.names = FALSE)

system(paste0("rm All_Train.txt"))
system(paste0("rm All_Tune.txt"))
system(paste0("rm All_Validation.txt"))