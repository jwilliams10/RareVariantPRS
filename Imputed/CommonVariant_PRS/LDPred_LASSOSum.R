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

time <- system.time({
  continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
  
  binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")

  NCORES <- 1
  
  trait <- as.numeric(commandArgs(TRUE)[1])
  
  trait <- c(continuous_traits,binary_traits)[trait]
  
  if(trait %in% continuous_traits){
    ## Continuous
    
    dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_",trait,".regenie"), sep="")
    colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
    dat$P <- 10^(-1*dat$LOG10P)
    
    sumstats <- dat[,c('CHROM', 'ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'P', 'N')]
    names(sumstats) <- c("chr", "rsid", "pos", "a1", "a0", "beta", "beta_se", "p", "n_eff")
    
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/all_chr_EUR_reference_Map.RData"))
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/all_chr_EUR_reference_Corr.RData"))
    
    sumstats <- sumstats[sumstats$rsid %in% map$rsid,]
    
    info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
    rownames(info_snp) <- info_snp$rsid
    
    df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
    
    print("Made it")
    
    corr <- corr[info_snp$`_NUM_ID_`,info_snp$`_NUM_ID_`]
    corr <- as_SFBM(as(corr, "generalMatrix"))
    
    print("Made it")
    
    # Automatic model
    ldsc <- snp_ldsc2(corr, df_beta)
    h2_est <- ldsc[["h2"]]
    print(paste0('Complete data preparation'))
    
    ## LDpred2
    h2_seq <- seq(0.1,1.5,by = 0.1)
    p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
    params <- expand.grid(p = p_seq, h2 = signif(abs(h2_est) * h2_seq, 3), sparse = c(FALSE))
    
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = NCORES)
    rownames(beta_grid) <- info_snp$rsid
    beta_grid <- cbind(beta_grid, info_snp[,c('a0','a1','rsid')])
    colnames(beta_grid) <- c(paste0('e',1:nrow(params)), 'a0','a1','rsid')
    beta_grid[is.na(beta_grid)] <- 0
    beta_grid <- as.data.frame(beta_grid)
    
    ## LASSOsum2
    delta_path <- function (max=100, min=0.5, n=10){
      sqrt_max <- max^(1/3)
      sqrt_min <- min^(1/3)
      path <- numeric(n)
      for (i in 1:n) {
        path[n+1-i] = (sqrt_max - (sqrt_max-sqrt_min)/(n-1)  * (i-1) )^3;
      }
      return(path)
    }
    beta_LASSOsum2 <- snp_lassosum2(corr, df_beta, delta = delta_path(max=100,min=0.5,n=10),ncores = NCORES, maxiter=1000)
    params2 <- attr(beta_LASSOsum2, "grid_param")
    rownames(beta_LASSOsum2) <- info_snp$rsid
    beta_LASSOsum2 <- cbind(beta_LASSOsum2, info_snp[,c('a0','a1','rsid')])
    colnames(beta_LASSOsum2) <- c(paste0('e',1:nrow(params2)), 'a0','a1','rsid')
    beta_LASSOsum2[is.na(beta_LASSOsum2)] <- 0
    beta_LASSOsum2 <- as.data.frame(beta_LASSOsum2)
    
    rm(corr)
    print(paste0('Complete'))
    
    # -------- PRS:
    
    ## LDpred2 
    prs.file <- data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, REF = beta_grid$a1, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
    write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt"),col.names = T,row.names = F,quote=F)
    
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation"))
    
    ## LASSOsum2
    prs.file <- data.frame(SNP = beta_LASSOsum2$rsid, ALT = beta_LASSOsum2$a0, REF = beta_LASSOsum2$a1, BETA = beta_LASSOsum2[,1:(ncol(beta_LASSOsum2)-3)])
    write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt"),col.names = T,row.names = F,quote=F)
    
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation"))
    
    ################ LDPred
    
    prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train.log")))
    prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune.log")))
    prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation.log")))
    
    sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))
    
    pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
    pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")
    
    pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
    pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")
    
    pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
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
    
    write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)
    
    ##### Final Coefficients
    
    # all_betas <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt"), sep="")
    # colnames(all_betas) <- c("SNP","ALT","REF",paste0("LDpred2_SCORE",1:nrow(sets),"_SUM"))
    # system(paste("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt")))
    # 
    # dat <- read.delim(paste0("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/GWAS_SummaryStatistics/",trait,"_sumstats.",trait,".glm.linear"), header=FALSE, comment.char="#")
    # colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
    # dat <- dat[dat$TEST == "ADD",]
    # dat <- dat[,c("CHROM","ID","REF","POS","A1")]
    # colnames(dat) <- c("CHR","SNP","REF","BP","A1")
    # dat <- left_join(dat,all_betas)
    # dat[is.na(dat)] <- 0
    # write.csv(all_betas,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Final_Coefficients.csv"),row.names = FALSE)
    
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
    
    write.csv(ldpred2_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"),row.names = FALSE)
    
    ############################### LASSOsum2
    
    prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train.log")))
    prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune.log")))
    prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation.log")))
    
    ## Pull in Phenotypes/Covariates 
    pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
    pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")
    
    pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
    pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")
    
    pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
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
    
    write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)
    
    ##### Final Coefficients
    # all_betas <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt"), sep="")
    # colnames(all_betas) <- c("SNP","ALT","REF",paste0("LASSOsum2_SCORE",1:300,"_SUM"))
    # system(paste("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt")))
    # 
    # dat <- read.delim(paste0("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/GWAS_SummaryStatistics/",trait,"_sumstats.",trait,".glm.linear"), header=FALSE, comment.char="#")
    # colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
    # dat <- dat[dat$TEST == "ADD",]
    # dat <- dat[,c("CHROM","ID","REF","POS","A1")]
    # colnames(dat) <- c("CHR","SNP","REF","BP","A1")
    # dat <- left_join(dat,all_betas)
    # dat[is.na(dat)] <- 0
    # write.csv(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Final_Coefficients.csv"),row.names = FALSE)
    # 
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
    
    write.csv(lassosum2_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"),row.names = FALSE)
    
  }else{
    ## Binary
    
    if(trait %in% c("Breast","Prostate")){
      fill <- "bp"
      confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }else{
      fill <- "act"
      confounders <-  paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }
    
    dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_",fill,"_",trait,".regenie"), sep="")
    colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
    dat$P <- 10^(-1*dat$LOG10P)
    sumstats <- dat[,c('CHROM', 'ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'P', 'N')]
    names(sumstats) <- c("chr", "rsid", "pos", "a1", "a0", "beta", "beta_se", "p", "n_eff")
    
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/all_chr_EUR_reference_Map.RData"))
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/all_chr_EUR_reference_Corr.RData"))
    
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
    
    ## LASSOsum2
    delta_path <- function (max=100, min=0.5, n=10){
      sqrt_max <- max^(1/3)
      sqrt_min <- min^(1/3)
      path <- numeric(n)
      for (i in 1:n) {
        path[n+1-i] = (sqrt_max - (sqrt_max-sqrt_min)/(n-1)  * (i-1) )^3;
      }
      return(path)
    }
    beta_LASSOsum2 <- snp_lassosum2(corr, df_beta, delta = delta_path(max=100,min=0.5,n=10),ncores = NCORES, maxiter=1000)
    params2 <- attr(beta_LASSOsum2, "grid_param")
    rownames(beta_LASSOsum2) = info_snp$rsid
    beta_LASSOsum2 = cbind(beta_LASSOsum2, info_snp[,c('a0','a1','rsid')])
    colnames(beta_LASSOsum2) = c(paste0('e',1:nrow(params2)), 'a0','a1','rsid')
    beta_LASSOsum2[is.na(beta_LASSOsum2)] <- 0
    beta_LASSOsum2 = as.data.frame(beta_LASSOsum2)
    
    rm(corr0, corr)
    print(paste0('Complete'))
    
    # -------- PRS:
    
    ## LDpred2 
    prs.file <- data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, REF = beta_grid$a1, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
    write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt"),col.names = T,row.names = F,quote=F)
    
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation"))
    
    ## LASSOsum2
    prs.file <- data.frame(SNP = beta_LASSOsum2$rsid, ALT = beta_LASSOsum2$a0, REF = beta_LASSOsum2$a1, BETA = beta_LASSOsum2[,1:(ncol(beta_LASSOsum2)-3)])
    write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt"),col.names = T,row.names = F,quote=F)
    
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune"))
    system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt cols=+scoresums,-scoreavgs header no-mean-imputation --score-col-nums 4-",ncol(prs.file)," --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation"))
    
    ################
    
    prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_train.log")))
    prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_tune.log")))
    prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_prs_validation.log")))
    
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
    
    for(k in 1:nrow(sets)){
      d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("SCORE",k,"_SUM"))]
      # d[,paste0("SCORE",k,"_SUM")] <- -1*d[,paste0("SCORE",k,"_SUM")]
      
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
    
    write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)
    
    ##### Final Coefficients
    # all_betas <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt"), sep="")
    # colnames(all_betas) <- c("SNP","ALT","REF",paste0("LDpred2_SCORE",1:nrow(sets),"_SUM"))
    # system(paste("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_LDpred2.txt")))
    # 
    # dat <- read.csv(paste0("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/GWAS_SummaryStatistics/regenie_step2_",fill,"_",trait,".regenie"), sep="")
    # colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
    # dat <- dat[,c("CHROM","ID","REF","POS","ALT")]
    # colnames(dat) <- c("CHR","SNP","REF","BP","A1")
    # dat <- left_join(dat,all_betas)
    # dat[is.na(dat)] <- 0
    # write.csv(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"_Final_Coefficients.csv"),row.names = FALSE)
    # 
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    pheno_validation_raw <- pheno_vad
    pheno_validation_adjusted <- pheno_vad
    
    mod <- lm(as.formula(paste0(paste0("SCORE",idx,"_SUM"),"~pc1 + pc2 + pc3 + pc4 + pc5")),data = pheno_validation_adjusted)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
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
    
    pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")])
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
      return(c(result))
    }
    
    AUC_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    beta_lower_validation_raw_EUR <- beta_ci$basic[4]
    beta_upper_validation_raw_EUR <- beta_ci$basic[5]
    
    auc_validation_raw_EUR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_EUR[!is.na(pheno_validation_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_EUR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_EUR <- sd(boot_auc$t)
    auc_lower_validation_raw_EUR <- auc_ci$basic[4]
    auc_upper_validation_raw_EUR <- auc_ci$basic[5]
    
    beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    beta_lower_validation_raw_SAS <- beta_ci$basic[4]
    beta_upper_validation_raw_SAS <- beta_ci$basic[5]
    
    auc_validation_raw_SAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_SAS[!is.na(pheno_validation_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_SAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_SAS <- sd(boot_auc$t)
    auc_lower_validation_raw_SAS <- auc_ci$basic[4]
    auc_upper_validation_raw_SAS <- auc_ci$basic[5]
    
    beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    beta_lower_validation_raw_AMR <- beta_ci$basic[4]
    beta_upper_validation_raw_AMR <- beta_ci$basic[5]
    
    auc_validation_raw_AMR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_AMR[!is.na(pheno_validation_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_AMR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_AMR <- sd(boot_auc$t)
    auc_lower_validation_raw_AMR <- auc_ci$basic[4]
    auc_upper_validation_raw_AMR <- auc_ci$basic[5]
    
    beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    beta_lower_validation_raw_AFR <- beta_ci$basic[4]
    beta_upper_validation_raw_AFR <- beta_ci$basic[5]
    
    auc_validation_raw_AFR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_AFR[!is.na(pheno_validation_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_AFR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_AFR <- sd(boot_auc$t)
    auc_lower_validation_raw_AFR <- auc_ci$basic[4]
    auc_upper_validation_raw_AFR <- auc_ci$basic[5]
    
    beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_EAS <- sd(boot_beta$t)
    beta_lower_validation_raw_EAS <- beta_ci$basic[4]
    beta_upper_validation_raw_EAS <- beta_ci$basic[5]
    
    auc_validation_raw_EAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_EAS[!is.na(pheno_validation_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_EAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_EAS <- sd(boot_auc$t)
    auc_lower_validation_raw_EAS <- auc_ci$basic[4]
    auc_upper_validation_raw_EAS <- auc_ci$basic[5]
    
    beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
    beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]
    
    auc_validation_adjusted_EUR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_EUR[!is.na(pheno_validation_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_EUR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_EUR <- sd(boot_auc$t)
    auc_lower_validation_adjusted_EUR <- auc_ci$basic[4]
    auc_upper_validation_adjusted_EUR <- auc_ci$basic[5]
    
    beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
    beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]
    
    auc_validation_adjusted_SAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_SAS[!is.na(pheno_validation_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_SAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_SAS <- sd(boot_auc$t)
    auc_lower_validation_adjusted_SAS <- auc_ci$basic[4]
    auc_upper_validation_adjusted_SAS <- auc_ci$basic[5]
    
    beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
    beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]
    
    auc_validation_adjusted_AMR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_AMR[!is.na(pheno_validation_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_AMR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_AMR <- sd(boot_auc$t)
    auc_lower_validation_adjusted_AMR <- auc_ci$basic[4]
    auc_upper_validation_adjusted_AMR <- auc_ci$basic[5]
    
    beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
    beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]
    
    auc_validation_adjusted_AFR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_AFR[!is.na(pheno_validation_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_AFR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_AFR <- sd(boot_auc$t)
    auc_lower_validation_adjusted_AFR <- auc_ci$basic[4]
    auc_upper_validation_adjusted_AFR <- auc_ci$basic[5]
    
    beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
    beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
    beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]
    
    auc_validation_adjusted_EAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_EAS[!is.na(pheno_validation_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_EAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_EAS <- sd(boot_auc$t)
    auc_lower_validation_adjusted_EAS <- auc_ci$basic[4]
    auc_upper_validation_adjusted_EAS <- auc_ci$basic[5]
    
    ldpred2_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                  beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                  beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                                  beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR,beta_lower_validation_raw_EAS), 
                                  beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR,beta_upper_validation_raw_EAS), 
                                  AUC_raw = c(auc_validation_raw_EUR,auc_validation_raw_SAS,auc_validation_raw_AMR,auc_validation_raw_AFR,auc_validation_raw_EAS),
                                  AUC_se_raw = c(auc_se_validation_raw_EUR,auc_se_validation_raw_SAS,auc_se_validation_raw_AMR,auc_se_validation_raw_AFR,auc_se_validation_raw_EAS),
                                  AUC_lower_raw = c(auc_lower_validation_raw_EUR,auc_lower_validation_raw_SAS,auc_lower_validation_raw_AMR,auc_lower_validation_raw_AFR,auc_lower_validation_raw_EAS),
                                  AUC_upper_raw = c(auc_upper_validation_raw_EUR,auc_upper_validation_raw_SAS,auc_upper_validation_raw_AMR,auc_upper_validation_raw_AFR,auc_upper_validation_raw_EAS),
                                  beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                  beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                                  beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR,beta_lower_validation_adjusted_EAS), 
                                  beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR,beta_upper_validation_adjusted_EAS), 
                                  AUC_adjusted = c(auc_validation_adjusted_EUR,auc_validation_adjusted_SAS,auc_validation_adjusted_AMR,auc_validation_adjusted_AFR,auc_validation_adjusted_EAS),
                                  AUC_se_adjusted = c(auc_se_validation_adjusted_EUR,auc_se_validation_adjusted_SAS,auc_se_validation_adjusted_AMR,auc_se_validation_adjusted_AFR,auc_se_validation_adjusted_EAS),
                                  AUC_lower_adjusted = c(auc_lower_validation_adjusted_EUR,auc_lower_validation_adjusted_SAS,auc_lower_validation_adjusted_AMR,auc_lower_validation_adjusted_AFR,auc_lower_validation_adjusted_EAS),
                                  AUC_upper_adjusted = c(auc_upper_validation_adjusted_EUR,auc_upper_validation_adjusted_SAS,auc_upper_validation_adjusted_AMR,auc_upper_validation_adjusted_AFR,auc_upper_validation_adjusted_EAS))
    
    write.csv(ldpred2_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"),row.names = FALSE)
    
    
    ########### LASSOsum2
    prs_mat_train <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_train.log")))
    prs_mat_tune <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_tune.log")))
    prs_mat_validation <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation.sscore"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_prs_validation.log")))
    
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
      # d[,paste0("SCORE",k,"_SUM")] <- -1*d[,paste0("SCORE",k,"_SUM")]
      
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
    
    write.table(best_prs_train,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2_train_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_tune,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2_tune_prs_best.txt"),sep = "\t",row.names = FALSE)
    write.table(best_prs_validation,file=paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2_validation_prs_best.txt"),sep = "\t",row.names = FALSE)
    
    ##### Final Coefficients
    # all_betas <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt"), sep="")
    # colnames(all_betas) <- c("SNP","ALT","REF",paste0("LASSOsum2_SCORE",1:300,"_SUM"))
    # system(paste("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_LASSOsum2.txt")))
    # 
    # dat <- read.csv(paste0("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/GWAS_SummaryStatistics/regenie_step2_",fill,"_",trait,".regenie"), sep="")
    # colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
    # dat <- dat[,c("CHROM","ID","REF","POS","ALT")]
    # colnames(dat) <- c("CHR","SNP","REF","BP","A1")
    # dat <- left_join(dat,all_betas)
    # dat[is.na(dat)] <- 0
    # write.csv(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"_Final_Coefficients.csv"),row.names = FALSE)
    
    
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
    pheno_validation_raw_SAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_validation_raw_AMR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    pheno_validation_raw_AFR <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_validation_raw_EAS <- pheno_validation_raw[pheno_validation_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    pheno_validation_adjusted_EUR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    pheno_validation_adjusted_SAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    pheno_validation_adjusted_AMR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    pheno_validation_adjusted_AFR <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    pheno_validation_adjusted_EAS <- pheno_validation_adjusted[pheno_validation_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EUR[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_SAS[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AMR[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_AFR[,paste0("SCORE",idx,"_SUM")])
    pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")] <- scale(pheno_validation_raw_EAS[,paste0("SCORE",idx,"_SUM")])
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
      return(c(result))
    }
    
    AUC_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    beta_lower_validation_raw_EUR <- beta_ci$basic[4]
    beta_upper_validation_raw_EUR <- beta_ci$basic[5]
    
    auc_validation_raw_EUR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_EUR[!is.na(pheno_validation_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_EUR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_EUR <- sd(boot_auc$t)
    auc_lower_validation_raw_EUR <- auc_ci$basic[4]
    auc_upper_validation_raw_EUR <- auc_ci$basic[5]
    
    beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    beta_lower_validation_raw_SAS <- beta_ci$basic[4]
    beta_upper_validation_raw_SAS <- beta_ci$basic[5]
    
    auc_validation_raw_SAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_SAS[!is.na(pheno_validation_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_SAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_SAS <- sd(boot_auc$t)
    auc_lower_validation_raw_SAS <- auc_ci$basic[4]
    auc_upper_validation_raw_SAS <- auc_ci$basic[5]
    
    beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    beta_lower_validation_raw_AMR <- beta_ci$basic[4]
    beta_upper_validation_raw_AMR <- beta_ci$basic[5]
    
    auc_validation_raw_AMR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_AMR[!is.na(pheno_validation_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_AMR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_AMR <- sd(boot_auc$t)
    auc_lower_validation_raw_AMR <- auc_ci$basic[4]
    auc_upper_validation_raw_AMR <- auc_ci$basic[5]
    
    beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    beta_lower_validation_raw_AFR <- beta_ci$basic[4]
    beta_upper_validation_raw_AFR <- beta_ci$basic[5]
    
    auc_validation_raw_AFR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_AFR[!is.na(pheno_validation_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_AFR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_AFR <- sd(boot_auc$t)
    auc_lower_validation_raw_AFR <- auc_ci$basic[4]
    auc_upper_validation_raw_AFR <- auc_ci$basic[5]
    
    beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_raw_EAS <- sd(boot_beta$t)
    beta_lower_validation_raw_EAS <- beta_ci$basic[4]
    beta_upper_validation_raw_EAS <- beta_ci$basic[5]
    
    auc_validation_raw_EAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_raw_EAS[!is.na(pheno_validation_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_raw_EAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_raw_EAS <- sd(boot_auc$t)
    auc_lower_validation_raw_EAS <- auc_ci$basic[4]
    auc_upper_validation_raw_EAS <- auc_ci$basic[5]
    
    beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
    beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]
    
    auc_validation_adjusted_EUR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_EUR[!is.na(pheno_validation_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_EUR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_EUR <- sd(boot_auc$t)
    auc_lower_validation_adjusted_EUR <- auc_ci$basic[4]
    auc_upper_validation_adjusted_EUR <- auc_ci$basic[5]
    
    beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
    beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]
    
    auc_validation_adjusted_SAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_SAS[!is.na(pheno_validation_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_SAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_SAS <- sd(boot_auc$t)
    auc_lower_validation_adjusted_SAS <- auc_ci$basic[4]
    auc_upper_validation_adjusted_SAS <- auc_ci$basic[5]
    
    beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
    beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]
    
    auc_validation_adjusted_AMR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_AMR[!is.na(pheno_validation_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_AMR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_AMR <- sd(boot_auc$t)
    auc_lower_validation_adjusted_AMR <- auc_ci$basic[4]
    auc_upper_validation_adjusted_AMR <- auc_ci$basic[5]
    
    beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
    beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]
    
    auc_validation_adjusted_AFR <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_AFR[!is.na(pheno_validation_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_AFR, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_AFR <- sd(boot_auc$t)
    auc_lower_validation_adjusted_AFR <- auc_ci$basic[4]
    auc_upper_validation_adjusted_AFR <- auc_ci$basic[5]
    
    beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~",paste0("SCORE",idx,"_SUM"),"+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 1000)
    beta_ci <- boot.ci(boot_beta, type = "basic")
    beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
    beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
    beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]
    
    auc_validation_adjusted_EAS <- roc.binary(status = trait,variable = paste0("SCORE",idx,"_SUM"),confounders = as.formula(confounders),data = pheno_validation_adjusted_EAS[!is.na(pheno_validation_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_auc <- boot(data = pheno_validation_adjusted_EAS, statistic = AUC_Boot, R = 1000)
    auc_ci <- boot.ci(boot_auc, type = "basic")
    auc_se_validation_adjusted_EAS <- sd(boot_auc$t)
    auc_lower_validation_adjusted_EAS <- auc_ci$basic[4]
    auc_upper_validation_adjusted_EAS <- auc_ci$basic[5]
    
    lassosum2_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                                    beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                                    beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                                    beta_lower_raw = c(beta_lower_validation_raw_EUR,beta_lower_validation_raw_SAS,beta_lower_validation_raw_AMR,beta_lower_validation_raw_AFR,beta_lower_validation_raw_EAS), 
                                    beta_upper_raw = c(beta_upper_validation_raw_EUR,beta_upper_validation_raw_SAS,beta_upper_validation_raw_AMR,beta_upper_validation_raw_AFR,beta_upper_validation_raw_EAS), 
                                    AUC_raw = c(auc_validation_raw_EUR,auc_validation_raw_SAS,auc_validation_raw_AMR,auc_validation_raw_AFR,auc_validation_raw_EAS),
                                    AUC_se_raw = c(auc_se_validation_raw_EUR,auc_se_validation_raw_SAS,auc_se_validation_raw_AMR,auc_se_validation_raw_AFR,auc_se_validation_raw_EAS),
                                    AUC_lower_raw = c(auc_lower_validation_raw_EUR,auc_lower_validation_raw_SAS,auc_lower_validation_raw_AMR,auc_lower_validation_raw_AFR,auc_lower_validation_raw_EAS),
                                    AUC_upper_raw = c(auc_upper_validation_raw_EUR,auc_upper_validation_raw_SAS,auc_upper_validation_raw_AMR,auc_upper_validation_raw_AFR,auc_upper_validation_raw_EAS),
                                    beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                                    beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                                    beta_lower_adjusted = c(beta_lower_validation_adjusted_EUR,beta_lower_validation_adjusted_SAS,beta_lower_validation_adjusted_AMR,beta_lower_validation_adjusted_AFR,beta_lower_validation_adjusted_EAS), 
                                    beta_upper_adjusted = c(beta_upper_validation_adjusted_EUR,beta_upper_validation_adjusted_SAS,beta_upper_validation_adjusted_AMR,beta_upper_validation_adjusted_AFR,beta_upper_validation_adjusted_EAS), 
                                    AUC_adjusted = c(auc_validation_adjusted_EUR,auc_validation_adjusted_SAS,auc_validation_adjusted_AMR,auc_validation_adjusted_AFR,auc_validation_adjusted_EAS),
                                    AUC_se_adjusted = c(auc_se_validation_adjusted_EUR,auc_se_validation_adjusted_SAS,auc_se_validation_adjusted_AMR,auc_se_validation_adjusted_AFR,auc_se_validation_adjusted_EAS),
                                    AUC_lower_adjusted = c(auc_lower_validation_adjusted_EUR,auc_lower_validation_adjusted_SAS,auc_lower_validation_adjusted_AMR,auc_lower_validation_adjusted_AFR,auc_lower_validation_adjusted_EAS),
                                    AUC_upper_adjusted = c(auc_upper_validation_adjusted_EUR,auc_upper_validation_adjusted_SAS,auc_upper_validation_adjusted_AMR,auc_upper_validation_adjusted_AFR,auc_upper_validation_adjusted_EAS))
    
    write.csv(lassosum2_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"),row.names = FALSE)
    
  }
})[3]

save(time, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/",trait,"_LDpred2_Lassosum2_Time.RData"))

