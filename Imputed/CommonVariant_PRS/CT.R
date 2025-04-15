rm(list = ls())
library(data.table)
library(dplyr)
library(RISCA)
library(boot)

time <- system.time({
  continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
  
  binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")
  
  trait <- as.numeric(commandArgs(TRUE)[1])
  
  trait <- c(continuous_traits,binary_traits)[trait]
  
  if(trait %in% continuous_traits){
    
    dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_",trait,".regenie"), sep="")
    colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
    dat$P <- 10^(-1*dat$LOG10P)
    
    dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
    colnames(dat) <- c("CHR","SNP","REF","BP","ALT","BETA","P")
    
    write.table(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/",trait,"_assoc.txt"),col.names = T,row.names = F,quote=F)
    
    pthr <- 1
    r2thr <- 0.1
    kbpthr <- 500
    # Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
    system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference --clump /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/",trait,"_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ","/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump"))
    
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.nosex"))
    
    ################################################
    
    ## This is the beginning of the thresholding step of clumping and threshold
    
    #p-value thresholds
    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
    LD <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.clumped")))
    # grab the index SNP for each clump
    clump.snp <- LD[,3,drop=F]
    # join against the summary data, this is now a dataset with n = number of index SNPs
    prs.all <- left_join(clump.snp,dat)
    
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.clumped")))
    
    n_pthres <- length(pthres)
    
    ## Write Coefficients of index SNPs to use later
    prs.file <- prs.all[,c("SNP","ALT","BETA")]
    write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff"),col.names = T,row.names = F,quote=F)
    
    # Write p-values to file to use later
    p.value.file <- prs.all[,c("SNP","P")]
    write.table(p.value.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value"),col.names = T,row.names = F,quote=F)
    
    beta_final <- inner_join(prs.file,p.value.file)
    
    ## Write a file with the p-value thresholds 
    q_range <- data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
    temp <- 1
    for(k in 1:length(pthres)){
      q_range[temp,1] <- paste0("p_value_",k)
      q_range[temp,3] <- pthres[k]
      temp <- temp + 1
    }
    q_range <- q_range[1:(temp-1),]
    write.table(q_range,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file"),row.names = F,col.names = F,quote=F)
    
    #PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
    #PRS = G*beta
    # Literally just multiply the score file (weighted or unweighted coefficients) by the G matrix, q-score-range is only for C + T, for LD pred score file would be weight coefficients 
    
    system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value header --threads 2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train"))
    system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value header --threads 2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune"))
    system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value header --threads 2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation"))
    #########################################################################
    
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value")))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff")))
    
    ### Merge all the outputted files from the previous command into one large data file
    prs_list <- list()
    temp <- 1
    for(k in 1:length(pthres)){
      
      prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train.p_value_",k,".sscore"))
      system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train.p_value_",k,".sscore")))
      # times (2*number of SNPs)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      
      colnames(prs_list[[temp]]) <- paste0("p_value_",k)
      temp <- temp + 1
      
    }
    prs_mat_train <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
    colnames(prs_mat_train)[2] <- "IID"
    
    write.table(prs_mat_train,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_train.txt"),row.names = F)
    
    prs_list <- list()
    temp <- 1
    for(k in 1:length(pthres)){
      
      prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune.p_value_",k,".sscore"))
      system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune.p_value_",k,".sscore")))
      # times (2*number of SNPs)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      
      colnames(prs_list[[temp]]) <- paste0("p_value_",k)
      temp <- temp + 1
      
    }
    prs_mat_tune <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
    colnames(prs_mat_tune)[2] <- "IID"
    
    write.table(prs_mat_tune,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_tune.txt"),row.names = F)
    
    ### Merge all the outputted files from the previous command into one large data file
    prs_list <- list()
    temp <- 1
    for(k in 1:length(pthres)){
      
      prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation.p_value_",k,".sscore"))
      system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation.p_value_",k,".sscore")))
      # times (2*number of SNPs)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      
      colnames(prs_list[[temp]]) <- paste0("p_value_",k)
      temp <- temp + 1
      
    }
    prs_mat_validation <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
    colnames(prs_mat_validation)[2] <- "IID"
    
    write.table(prs_mat_validation,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_validation.txt"),row.names = F)
    
    
    ## Pull in Phenotypes/Covariates 
    pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
    pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")
    
    pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
    pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")
    
    pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
    pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")
    
    #calculate R2 for each of the tuning dataset
    # This is done by regressing the residuals of the model with all covariates against the prs
    
    r2_tun_vec <- rep(0,length(pthres))
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tuning)
    for(k in 1:length(pthres)){
      prs <- pheno_tuning[!is.na(pheno_tuning[,trait]),paste0("p_value_",k)]
      model.prs <- lm(model.null$residual~prs,data=pheno_tuning)
      r2_tun_vec[k] <- summary(model.prs)$r.square
    }
    
    #find best p-value threshold
    idx <- which.max(r2_tun_vec)
    #write down best prs in the prs folder/ save it and store it
    prs_train_max <- pheno_train[,c("IID","FID",paste0("p_value_",idx))]
    colnames(prs_train_max) <- c("IID","FID","prs")
    write.table(prs_train_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train_best.txt"),row.names = F)
    
    prs_tune_max <- pheno_tuning[,c("IID","FID",paste0("p_value_",idx))]
    colnames(prs_tune_max) <- c("IID","FID","prs")
    write.table(prs_tune_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune_best.txt"),row.names = F)
    
    prs_vad_max <- pheno_vad[,c("IID","FID",paste0("p_value_",idx))]
    colnames(prs_vad_max) <- c("IID","FID","prs")
    write.table(prs_vad_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation_best.txt"),row.names = F)
    
    ##### Final Coefficients
    # dat <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/",trait,"_sumstats.",trait,".glm.linear"), header=FALSE, comment.char="#")
    # colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
    # dat <- dat[dat$TEST == "ADD",]
    # dat <- dat[,c("CHROM","ID","REF","POS","ALT","P","BETA")]
    # colnames(dat) <- c("CHR","SNP","REF","BP","ALT","P","BETA")
    # 
    # dat$CT_p_value_1 <- 0
    # dat$CT_p_value_1[dat$P <= pthres[1]] <- dat$BETA[dat$P <= pthres[1]]
    # dat$CT_p_value_1[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_2 <- 0
    # dat$CT_p_value_2[dat$P <= pthres[2]] <- dat$BETA[dat$P <= pthres[2]]
    # dat$CT_p_value_2[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_3 <- 0
    # dat$CT_p_value_3[dat$P <= pthres[3]] <- dat$BETA[dat$P <= pthres[3]]
    # dat$CT_p_value_3[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_4 <- 0
    # dat$CT_p_value_4[dat$P <= pthres[4]] <- dat$BETA[dat$P <= pthres[4]]
    # dat$CT_p_value_4[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_5 <- 0
    # dat$CT_p_value_5[dat$P <= pthres[5]] <- dat$BETA[dat$P <= pthres[5]]
    # dat$CT_p_value_5[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_6 <- 0
    # dat$CT_p_value_6[dat$P <= pthres[6]] <- dat$BETA[dat$P <= pthres[6]]
    # dat$CT_p_value_6[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_7 <- 0
    # dat$CT_p_value_7[dat$P <= pthres[7]] <- dat$BETA[dat$P <= pthres[7]]
    # dat$CT_p_value_7[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_8 <- 0
    # dat$CT_p_value_8[dat$P <= pthres[8]] <- dat$BETA[dat$P <= pthres[8]]
    # dat$CT_p_value_8[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_9 <- 0
    # dat$CT_p_value_9[dat$P <= pthres[9]] <- dat$BETA[dat$P <= pthres[9]]
    # dat$CT_p_value_9[!(dat$SNP %in% beta_final$SNP)] <- 0
    # 
    # dat <- dat[,c("CHR","SNP","REF","BP","ALT","P",paste0("CT_p_value_",1:9))]
    # write.csv(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Final_Coefficients.csv"),row.names = FALSE)
    write.csv(beta_final,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_CT_SNP_List.csv"),row.names = FALSE)
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    pheno_vad$y_validation <- NA
    pheno_vad$y_validation[!is.na(pheno_vad[,trait])] <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_vad)$residual
    
    pheno_validation_raw <- pheno_vad
    pheno_validation_adjusted <- pheno_vad
    
    mod <- lm(as.formula(paste0(paste0("p_value_",idx),"~pc1 + pc2 + pc3 + pc4 + pc5")),data = pheno_validation_adjusted)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(sqrt(y_hat)) == 0){
      pheno_validation_adjusted[,paste0("p_value_",idx)] <- 0
    }else{
      pheno_validation_adjusted[,paste0("p_value_",idx)] <- R/sqrt(y_hat)
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
    
    pheno_validation_raw_EUR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EUR[,paste0("p_value_",idx)])
    pheno_validation_raw_SAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_SAS[,paste0("p_value_",idx)])
    pheno_validation_raw_AMR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AMR[,paste0("p_value_",idx)])
    pheno_validation_raw_AFR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AFR[,paste0("p_value_",idx)])
    pheno_validation_raw_EAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EAS[,paste0("p_value_",idx)])
    
    pheno_validation_adjusted_EUR$y_validation <- scale(pheno_validation_adjusted_EUR$y_validation)
    pheno_validation_adjusted_SAS$y_validation <- scale(pheno_validation_adjusted_SAS$y_validation)
    pheno_validation_adjusted_AMR$y_validation <- scale(pheno_validation_adjusted_AMR$y_validation)
    pheno_validation_adjusted_AFR$y_validation <- scale(pheno_validation_adjusted_AFR$y_validation)
    pheno_validation_adjusted_EAS$y_validation <- scale(pheno_validation_adjusted_EAS$y_validation)
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = boot_data))[2]
      return(c(result))
    }
    
    R2_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = boot_data))$r.squared
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_EUR))[2]
    boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 10000)
    beta_raw_EUR_boot <- boot_beta$t
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    
    R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_EUR))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_EUR, statistic = R2_Boot, R = 10000)
    R2_raw_EUR_boot <- boot_R2$t
    R2_se_validation_raw_EUR <- sd(boot_R2$t)
    
    beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_SAS))[2]
    boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 10000)
    beta_raw_SAS_boot <- boot_beta$t
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    
    R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_SAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_SAS, statistic = R2_Boot, R = 10000)
    R2_raw_SAS_boot <- boot_R2$t
    R2_se_validation_raw_SAS <- sd(boot_R2$t)
    
    beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AMR))[2]
    boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 10000)
    beta_raw_AMR_boot <- boot_beta$t
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    
    R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AMR))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_AMR, statistic = R2_Boot, R = 10000)
    R2_raw_AMR_boot <- boot_R2$t
    R2_se_validation_raw_AMR <- sd(boot_R2$t)
    
    beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AFR))[2]
    boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 10000)
    beta_raw_AFR_boot <- boot_beta$t
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    
    R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_AFR))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_AFR, statistic = R2_Boot, R = 10000)
    R2_raw_AFR_boot <- boot_R2$t
    R2_se_validation_raw_AFR <- sd(boot_R2$t)
    
    beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_EAS))[2]
    boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 10000)
    beta_raw_EAS_boot <- boot_beta$t
    beta_se_validation_raw_EAS <- sd(boot_beta$t)
    
    R2_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_raw_EAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_raw_EAS, statistic = R2_Boot, R = 10000)
    R2_raw_EAS_boot <- boot_R2$t
    R2_se_validation_raw_EAS <- sd(boot_R2$t)
    
    beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_EUR))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EUR_boot <- boot_beta$t
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_EUR))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_EUR, statistic = R2_Boot, R = 10000)
    R2_adjusted_EUR_boot <- boot_R2$t
    R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
    
    beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_SAS))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_SAS_boot <- boot_beta$t
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_SAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_SAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_SAS_boot <- boot_R2$t
    R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
    
    beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AMR))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AMR_boot <- boot_beta$t
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AMR))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_AMR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AMR_boot <- boot_R2$t
    R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
    
    beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AFR))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AFR_boot <- boot_beta$t
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_AFR))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_AFR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AFR_boot <- boot_R2$t
    R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
    
    beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_EAS))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EAS_boot <- boot_beta$t
    beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~",paste0("p_value_",idx))),data = pheno_validation_adjusted_EAS))$r.squared
    boot_R2 <- boot(data = pheno_validation_adjusted_EAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_EAS_boot <- boot_R2$t
    R2_se_validation_adjusted_EAS <- sd(boot_R2$t)
    
    CT_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                             beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                             R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR,R2_validation_raw_EAS),
                             R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR,R2_se_validation_raw_EAS),
                             beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                             beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                             R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR,R2_validation_adjusted_EAS),
                             R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR,R2_se_validation_adjusted_EAS))
    
    CT_Boot_Results <- data.frame(trait = trait,beta_raw_EUR_boot,R2_raw_EUR_boot,beta_raw_SAS_boot,R2_raw_SAS_boot,
                                  beta_raw_AMR_boot,R2_raw_AMR_boot,beta_raw_AFR_boot,R2_raw_AFR_boot,
                                  beta_raw_EAS_boot,R2_raw_EAS_boot,beta_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                                  beta_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                                  beta_adjusted_AFR_boot,R2_adjusted_AFR_boot,beta_adjusted_EAS_boot,R2_adjusted_EAS_boot)
    
    write.csv(CT_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"),row.names = FALSE)
    write.csv(CT_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Bootstraps.csv"),row.names = FALSE)  
  }else{
    
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
    
    dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
    colnames(dat) <- c("CHR","SNP","REF","BP","ALT","BETA","P")
    
    write.table(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/",trait,"_assoc.txt"),col.names = T,row.names = F,quote=F)
    
    pthr <- 1
    r2thr <- 0.1
    kbpthr <- 500
    # Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
    system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference --clump /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/",trait,"_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ","/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump"))
    
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.nosex"))
    
    ################################################
    
    ## This is the beginning of the thresholding step of clumping and threshold
    
    #p-value thresholds
    pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
    LD <- as.data.frame(fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.clumped")))
    # grab the index SNP for each clump
    clump.snp <- LD[,3,drop=F]
    # join against the summary data, this is now a dataset with n = number of index SNPs
    prs.all <- left_join(clump.snp,dat)
    
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_LDclump.clumped")))
    
    n_pthres <- length(pthres)
    
    ## Write Coefficients of index SNPs to use later
    prs.file <- prs.all[,c("SNP","ALT","BETA")]
    write.table(prs.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff"),col.names = T,row.names = F,quote=F)
    
    # Write p-values to file to use later
    p.value.file <- prs.all[,c("SNP","P")]
    write.table(p.value.file,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value"),col.names = T,row.names = F,quote=F)
    
    beta_final <- inner_join(prs.file,p.value.file)
    
    ## Write a file with the p-value thresholds 
    q_range <- data.frame(rep("p_value",n_pthres),rep(0,n_pthres),rep(0.5,n_pthres))
    temp <- 1
    for(k in 1:length(pthres)){
      q_range[temp,1] <- paste0("p_value_",k)
      q_range[temp,3] <- pthres[k]
      temp <- temp + 1
    }
    q_range <- q_range[1:(temp-1),]
    write.table(q_range,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file"),row.names = F,col.names = F,quote=F)
    
    #PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
    #PRS = G*beta
    # Literally just multiply the score file (weighted or unweighted coefficients) by the G matrix, q-score-range is only for C + T, for LD pred score file would be weight coefficients 
    
    system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value header --threads 2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train"))
    system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value header --threads 2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune"))
    system(paste0("/data/williamsjacr/software/plink2 --q-score-range /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value header --threads 2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation"))
    #########################################################################
    
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train.log"))
    system(paste0("rm /data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"q_range_file"))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_p_value")))
    system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_coeff")))
    
    ### Merge all the outputted files from the previous command into one large data file
    prs_list <- list()
    temp <- 1
    for(k in 1:length(pthres)){
      
      prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train.p_value_",k,".sscore"))
      system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train.p_value_",k,".sscore")))
      # times (2*number of SNPs)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      
      colnames(prs_list[[temp]]) <- paste0("p_value_",k)
      temp <- temp + 1
      
    }
    prs_mat_train <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
    colnames(prs_mat_train)[2] <- "IID"
    
    write.table(prs_mat_train,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_train.txt"),row.names = F)
    
    prs_list <- list()
    temp <- 1
    for(k in 1:length(pthres)){
      
      prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune.p_value_",k,".sscore"))
      system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune.p_value_",k,".sscore")))
      # times (2*number of SNPs)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      
      colnames(prs_list[[temp]]) <- paste0("p_value_",k)
      temp <- temp + 1
      
    }
    prs_mat_tune <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
    colnames(prs_mat_tune)[2] <- "IID"
    
    write.table(prs_mat_tune,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_tune.txt"),row.names = F)
    
    ### Merge all the outputted files from the previous command into one large data file
    prs_list <- list()
    temp <- 1
    for(k in 1:length(pthres)){
      
      prs_temp <- fread(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation.p_value_",k,".sscore"))
      system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation.p_value_",k,".sscore")))
      # times (2*number of SNPs)
      prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
      
      colnames(prs_list[[temp]]) <- paste0("p_value_",k)
      temp <- temp + 1
      
    }
    prs_mat_validation <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
    colnames(prs_mat_validation)[2] <- "IID"
    
    write.table(prs_mat_validation,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_all_validation.txt"),row.names = F)
    
    
    ## Pull in Phenotypes/Covariates 
    pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
    pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")
    
    pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
    pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")
    
    pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
    pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")
    
    #calculate AUC for each of the tuning dataset
    # This is done by regressing the residuals of the model with all covariates against the prs
    AUC_tun_vec <- rep(0,length(pthres))
    
    for(k in 1:length(pthres)){
      d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",k))]
      
      roc_obj <- roc.binary(status = trait,
                            variable = paste0("p_value_",k),
                            confounders = as.formula(confounders),
                            data = d,
                            precision=seq(0.05,0.95, by=0.05))
      
      AUC_tun_vec[k] <- roc_obj$auc
    }
    
    #find best p-value threshold
    idx <- which.max(AUC_tun_vec)
    #write down best prs in the prs folder/ save it and store it
    prs_train_max <- pheno_train[,c("IID","FID",paste0("p_value_",idx))]
    colnames(prs_train_max) <- c("IID","FID","prs")
    write.table(prs_train_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_train_best.txt"),row.names = F)
    
    prs_tune_max <- pheno_tuning[,c("IID","FID",paste0("p_value_",idx))]
    colnames(prs_tune_max) <- c("IID","FID","prs")
    write.table(prs_tune_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_tune_best.txt"),row.names = F)
    
    prs_vad_max <- pheno_vad[,c("IID","FID",paste0("p_value_",idx))]
    colnames(prs_vad_max) <- c("IID","FID","prs")
    write.table(prs_vad_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_prs_validation_best.txt"),row.names = F)
    
    
    ##### Final Coefficients
    # dat$CT_p_value_1 <- 0
    # dat$CT_p_value_1[dat$P <= pthres[1]] <- dat$BETA[dat$P <= pthres[1]]
    # dat$CT_p_value_1[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_2 <- 0
    # dat$CT_p_value_2[dat$P <= pthres[2]] <- dat$BETA[dat$P <= pthres[2]]
    # dat$CT_p_value_2[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_3 <- 0
    # dat$CT_p_value_3[dat$P <= pthres[3]] <- dat$BETA[dat$P <= pthres[3]]
    # dat$CT_p_value_3[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_4 <- 0
    # dat$CT_p_value_4[dat$P <= pthres[4]] <- dat$BETA[dat$P <= pthres[4]]
    # dat$CT_p_value_4[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_5 <- 0
    # dat$CT_p_value_5[dat$P <= pthres[5]] <- dat$BETA[dat$P <= pthres[5]]
    # dat$CT_p_value_5[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_6 <- 0
    # dat$CT_p_value_6[dat$P <= pthres[6]] <- dat$BETA[dat$P <= pthres[6]]
    # dat$CT_p_value_6[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_7 <- 0
    # dat$CT_p_value_7[dat$P <= pthres[7]] <- dat$BETA[dat$P <= pthres[7]]
    # dat$CT_p_value_7[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_8 <- 0
    # dat$CT_p_value_8[dat$P <= pthres[8]] <- dat$BETA[dat$P <= pthres[8]]
    # dat$CT_p_value_8[!(dat$SNP %in% beta_final$SNP)] <- 0
    # dat$CT_p_value_9 <- 0
    # dat$CT_p_value_9[dat$P <= pthres[9]] <- dat$BETA[dat$P <= pthres[9]]
    # dat$CT_p_value_9[!(dat$SNP %in% beta_final$SNP)] <- 0
    # 
    # dat <- dat[,c("CHR","SNP","REF","BP","A1","P",paste0("CT_p_value_",1:9))]
    # 
    # write.csv(dat,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Final_Coefficients.csv"),row.names = FALSE)
    write.csv(beta_final,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_CT_SNP_List.csv"),row.names = FALSE)
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    pheno_validation_raw <- pheno_vad
    pheno_validation_adjusted <- pheno_vad
    
    mod <- lm(as.formula(paste0(paste0("p_value_",idx),"~pc1 + pc2 + pc3 + pc4 + pc5")),data = pheno_validation_adjusted)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,pheno_validation_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(sqrt(y_hat)) == 0){
      pheno_validation_adjusted[,paste0("p_value_",idx)] <- 0
    }else{
      pheno_validation_adjusted[,paste0("p_value_",idx)] <- R/sqrt(y_hat)
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
    
    pheno_validation_raw_EUR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EUR[,paste0("p_value_",idx)])
    pheno_validation_raw_SAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_SAS[,paste0("p_value_",idx)])
    pheno_validation_raw_AMR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AMR[,paste0("p_value_",idx)])
    pheno_validation_raw_AFR[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_AFR[,paste0("p_value_",idx)])
    pheno_validation_raw_EAS[,paste0("p_value_",idx)] <- scale(pheno_validation_raw_EAS[,paste0("p_value_",idx)])
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
      return(c(result))
    }
    
    AUC_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EUR, statistic = Beta_Boot, R = 10000)
    beta_raw_EUR_boot <- boot_beta$t
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    
    AUC_validation_raw_EUR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_EUR[!is.na(pheno_validation_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_EUR, statistic = AUC_Boot, R = 10000)
    AUC_raw_EUR_boot <- boot_AUC$t
    AUC_se_validation_raw_EUR <- sd(boot_AUC$t)
    
    beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_SAS, statistic = Beta_Boot, R = 10000)
    beta_raw_SAS_boot <- boot_beta$t
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    
    AUC_validation_raw_SAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_SAS[!is.na(pheno_validation_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_SAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_SAS_boot <- boot_AUC$t
    AUC_se_validation_raw_SAS <- sd(boot_AUC$t)
    
    beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AMR, statistic = Beta_Boot, R = 10000)
    beta_raw_AMR_boot <- boot_beta$t
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    
    AUC_validation_raw_AMR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_AMR[!is.na(pheno_validation_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_AMR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AMR_boot <- boot_AUC$t
    AUC_se_validation_raw_AMR <- sd(boot_AUC$t)
    
    beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_AFR, statistic = Beta_Boot, R = 10000)
    beta_raw_AFR_boot <- boot_beta$t
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    
    AUC_validation_raw_AFR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_AFR[!is.na(pheno_validation_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_AFR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AFR_boot <- boot_AUC$t
    AUC_se_validation_raw_AFR <- sd(boot_AUC$t)
    
    beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_raw_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_raw_EAS, statistic = Beta_Boot, R = 10000)
    beta_raw_EAS_boot <- boot_beta$t
    beta_se_validation_raw_EAS <- sd(boot_beta$t)
    
    AUC_validation_raw_EAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_raw_EAS[!is.na(pheno_validation_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_raw_EAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_EAS_boot <- boot_AUC$t
    AUC_se_validation_raw_EAS <- sd(boot_AUC$t)
    
    beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_EUR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EUR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EUR_boot <- boot_beta$t
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_EUR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_EUR[!is.na(pheno_validation_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_EUR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EUR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_EUR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_SAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_SAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_SAS_boot <- boot_beta$t
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    AUC_validation_adjusted_SAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_SAS[!is.na(pheno_validation_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_SAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_SAS_boot <- boot_AUC$t
    AUC_se_validation_adjusted_SAS <- sd(boot_AUC$t)
    
    beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_AMR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AMR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AMR_boot <- boot_beta$t
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_AMR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_AMR[!is.na(pheno_validation_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_AMR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AMR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_AMR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_AFR,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_AFR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AFR_boot <- boot_beta$t
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_AFR <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_AFR[!is.na(pheno_validation_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_AFR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AFR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_AFR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~","prs","+",gsub("~","",confounders))),data = pheno_validation_adjusted_EAS,family = binomial()))[2]
    boot_beta <- boot(data = pheno_validation_adjusted_EAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EAS_boot <- boot_beta$t
    beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
    
    AUC_validation_adjusted_EAS <- roc.binary(status = trait,variable = "prs",confounders = as.formula(confounders),data = pheno_validation_adjusted_EAS[!is.na(pheno_validation_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
    boot_AUC <- boot(data = pheno_validation_adjusted_EAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EAS_boot <- boot_AUC$t
    AUC_se_validation_adjusted_EAS <- sd(boot_AUC$t)
    
    CT_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                             beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                             AUC_raw = c(AUC_validation_raw_EUR,AUC_validation_raw_SAS,AUC_validation_raw_AMR,AUC_validation_raw_AFR,AUC_validation_raw_EAS),
                             AUC_se_raw = c(AUC_se_validation_raw_EUR,AUC_se_validation_raw_SAS,AUC_se_validation_raw_AMR,AUC_se_validation_raw_AFR,AUC_se_validation_raw_EAS),
                             beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                             beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                             AUC_adjusted = c(AUC_validation_adjusted_EUR,AUC_validation_adjusted_SAS,AUC_validation_adjusted_AMR,AUC_validation_adjusted_AFR,AUC_validation_adjusted_EAS),
                             AUC_se_adjusted = c(AUC_se_validation_adjusted_EUR,AUC_se_validation_adjusted_SAS,AUC_se_validation_adjusted_AMR,AUC_se_validation_adjusted_AFR,AUC_se_validation_adjusted_EAS))
    
    CT_Boot_Results <- data.frame(trait = trait,beta_raw_EUR_boot,AUC_raw_EUR_boot,beta_raw_SAS_boot,AUC_raw_SAS_boot,
                                  beta_raw_AMR_boot,AUC_raw_AMR_boot,beta_raw_AFR_boot,AUC_raw_AFR_boot,
                                  beta_raw_EAS_boot,AUC_raw_EAS_boot,beta_adjusted_EUR_boot,AUC_adjusted_EUR_boot,
                                  beta_adjusted_SAS_boot,AUC_adjusted_SAS_boot,beta_adjusted_AMR_boot,AUC_adjusted_AMR_boot,
                                  beta_adjusted_AFR_boot,AUC_adjusted_AFR_boot,beta_adjusted_EAS_boot,AUC_adjusted_EAS_boot)
    
    write.csv(CT_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"),row.names = FALSE) 
    write.csv(CT_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Bootstraps.csv"),row.names = FALSE) 
  }
})[3]

save(time,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"_Time.RData"))