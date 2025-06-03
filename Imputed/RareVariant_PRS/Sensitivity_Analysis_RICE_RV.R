rm(list = ls())
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(dplyr)
library(RISCA)
library(boot)
library(stringr)
library(caret)
library(ranger)
library(glmnet)

index <- as.numeric(commandArgs(TRUE)[1])

if(index < 12){
  thresholds <- 1e-5
}else if(index < 23){
  thresholds <- 1e-4
  index <- index - 11
}else if(index < 34){
  thresholds <- 1e-3
  index <- index - 22
}else if(index < 45){
  thresholds <- 1e-2
  index <- index - 33
}else{
  thresholds <- 1e-1
  index <- index - 44
}

source("~/RareVariantPRS/Simulation_Study/FullData/G_star/g_star_gene_centric_coding.R")

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")

binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")

trait <- c(continuous_traits,binary_traits)[index]

Train_PVals_All <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))

if(sum(Train_PVals_All$STAARB <= thresholds) == 0){
  RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR"), 
                               beta_raw = c(0,0,0,0), 
                               beta_se_raw = c(0,0,0,0), 
                               R2_raw = c(0,0,0,0),
                               R2_se_raw = c(0,0,0,0),
                               beta_adjusted = c(0,0,0,0), 
                               beta_se_adjusted = c(0,0,0,0), 
                               R2_adjusted = c(0,0,0,0),
                               R2_se_adjusted = c(0,0,0,0))
  
  RV_Boot_Results <- data.frame(trait = trait,c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),
                                c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),
                                c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),
                                c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),
                                c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0),c(0,0,0,0,0))
}else{
  Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= thresholds,]
  
  ## agds dir
  agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))
  
  ## Null Model
  obj_nullmodel_train <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/nullmodels_staar/",trait,"_Train_Null_Model.RData")))
  
  obj_nullmodel <- obj_nullmodel_train
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  pheno_valid <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  obj_nullmodel$id_include <- c(pheno_train$IID,pheno_tune$IID,pheno_valid$IID)
  
  ## Parameter
  QC_label <- "annotation/info/QC_label"
  geno_missing_imputation <- "mean"
  variant_type <- "SNV"	
  
  ## Annotation_dir
  Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
  ## Annotation channel
  Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/Annotation_name_catalog.Rdata"))
  
  G_star_gene_centric_coding <- list()
  
  for(i in 1:nrow(Train_PVals_All)){
    ## Chr
    chr <- Train_PVals_All$Chr[i]
    ## Gene name
    gene_name <- Train_PVals_All$Gene[i]
    ## Coding mask
    category <- Train_PVals_All$Category[i]
    
    ### gds file
    gds.path <- agds_dir[chr]
    genofile <- seqOpen(gds.path)
    
    G_star_gene_centric_coding[[i]] <- Gene_Centric_Coding_G_Star(chr=chr,gene_name=gene_name,category=category ,
                                                                  genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE) 
    seqClose(genofile) 
  } 
  
  G_star_gene_centric_coding <- do.call(cbind,G_star_gene_centric_coding)
  
  G_star_gene_centric_coding <- round(G_star_gene_centric_coding)
  
  col_remove <- apply(G_star_gene_centric_coding,2,function(x){sum(x != 0)}) > 10 & colSums(G_star_gene_centric_coding) > 10
  
  G_star_gene_centric_coding <- G_star_gene_centric_coding[,col_remove,drop = FALSE]
  
  Train_PVals_All <- Train_PVals_All[col_remove,,drop = FALSE]
  
  ids_gstar <- obj_nullmodel$id_include
  
  G_star_gene_centric_coding_train <- G_star_gene_centric_coding[ids_gstar %in% pheno_train$IID,,drop = FALSE]
  G_star_gene_centric_coding_tune <- G_star_gene_centric_coding[ids_gstar %in% pheno_tune$IID,,drop = FALSE]
  G_star_gene_centric_coding_vad <- G_star_gene_centric_coding[ids_gstar %in% pheno_valid$IID,,drop = FALSE]
  
  X_train <- data.frame(IID = ids_gstar[ids_gstar %in% pheno_train$IID],G_star_gene_centric_coding_train)
  colnames(X_train) <- c("IID",paste0("X",1:ncol(G_star_gene_centric_coding_train)))
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  pheno_train <- inner_join(pheno_train,X_train)
  X_train <- as.matrix(pheno_train[,paste0("X",1:ncol(G_star_gene_centric_coding_train)),drop = FALSE])
  
  X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% pheno_tune$IID],G_star_gene_centric_coding_tune)
  colnames(X_tune) <- c("IID",paste0("X",1:ncol(G_star_gene_centric_coding_tune)))
  pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  pheno_tune <- inner_join(pheno_tune,X_tune)
  X_tune <- as.matrix(pheno_tune[,paste0("X",1:ncol(G_star_gene_centric_coding_tune)),drop = FALSE])
  
  X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% pheno_valid$IID],G_star_gene_centric_coding_vad)
  colnames(X_valid) <- c("IID",paste0("X",1:ncol(G_star_gene_centric_coding_vad)))
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  pheno_validation <- inner_join(pheno_validation,X_valid)
  X_valid <- as.matrix(pheno_validation[,paste0("X",1:ncol(G_star_gene_centric_coding_vad)),drop = FALSE])
  
  Ensemble_Function_Continuous <- function(x,y){
    x <- as.matrix(x[!is.na(y),])
    y <- y[!is.na(y)]
    
    lasso_train <- glmnet(x,y,family = "gaussian",alpha = 1)
    ridge_train <- glmnet(x,y,family = "gaussian",alpha = 0)
    
    lasso_prs_tune <- predict(lasso_train,x)
    ridge_prs_tune <- predict(ridge_train,x)
    
    all <- cbind(lasso_prs_tune,ridge_prs_tune)
    
    R2_Vector <- vector()
    for(i in 1:ncol(all)){
      tmp <- data.frame(y = y, x_try = all[,i])
      R2_Vector[i] <- summary(lm(y~x_try,data = tmp))$r.square
    }
    
    coefficients_x <- coef(lm(y~.,data.frame(y = all[,which.max(R2_Vector)],x)))
    return(list(Coefficients = coefficients_x))
  }
  Ensemble_Function_Binary <- function(x,y){
    x <- as.matrix(x[!is.na(y),])
    y <- y[!is.na(y)]
    
    lasso_train <- glmnet(x,y,family = "binomial",alpha = 1)
    ridge_train <- glmnet(x,y,family = "binomial",alpha = 0)
    
    lasso_prs_tune <- predict(lasso_train,x)
    ridge_prs_tune <- predict(ridge_train,x)
    
    all <- cbind(lasso_prs_tune,ridge_prs_tune)
    
    AUC_Vector <- vector()
    for(i in 1:ncol(all)){
      tmp <- data.frame(y = y, x_try = all[,i])
      roc_obj <- roc.binary(status = "y",
                            variable = "x_try",
                            confounders = "~1",
                            data = tmp,
                            precision=seq(0.05,0.95, by=0.05))
      AUC_Vector[i] <- roc_obj$auc
    }
    
    coefficients_x <- coef(lm(y~.,data.frame(y = all[,which.max(AUC_Vector)],x)))
    return(list(Coefficients = coefficients_x))
  }
  Ensemble_Function <- function(x,y,family = c("continuous","binary")){
    if(family == "continuous"){
      return(Ensemble_Function_Continuous(x,y))
    }else{
      return(Ensemble_Function_Binary(x,y))
    }
  }
  
  if(trait %in% binary_traits){
    
    if(ncol(X_train) == 1){
      lm_train <- lm.fit(cbind(1,X_train[!is.na(pheno_train[,trait]),]),pheno_train[!is.na(pheno_train[,trait]),trait])
      lm_train$coefficients[is.na(lm_train$coefficients)] <- 0
      
      lm_prs_tune <- as.numeric(cbind(1,X_tune)%*%matrix(lm_train$coefficients,ncol = 1))
      
      lm_prs_vad <- as.numeric(cbind(1,X_valid)%*%matrix(lm_train$coefficients,ncol = 1))
      
      PRS_Tune <- as.matrix(lm_prs_tune,ncol = 1)
      PRS_Validation <- as.matrix(lm_prs_vad,ncol = 1)
    }else{
      
      lasso_train <- glmnet(X_train[!is.na(pheno_train[,trait]),],pheno_train[!is.na(pheno_train[,trait]),trait],family = "binomial",alpha = 1)
      ridge_train <- glmnet(X_train[!is.na(pheno_train[,trait]),],pheno_train[!is.na(pheno_train[,trait]),trait],family = "binomial",alpha = 0)
      lm_train <- glm(y~.,data = data.frame(y = pheno_train[!is.na(pheno_train[,trait]),trait], X_train[!is.na(pheno_train[,trait]),]),family = binomial())
      
      beta_matrix <- as.data.frame(cbind(as.matrix(lasso_train$beta),as.matrix(ridge_train$beta),lm_train$coefficients[-1]))
      colnames(beta_matrix) <- c(paste0("lasso_prs",1:ncol(as.matrix(lasso_train$beta))),paste0("ridge_prs",1:ncol(as.matrix(ridge_train$beta))),paste0("lm_prs",1))
      beta_matrix <- cbind(Train_PVals_All[,c(1:4)],beta_matrix)
      
      
      lasso_prs_tune <- predict(lasso_train,X_tune,type = "link")
      ridge_prs_tune <- predict(ridge_train,X_tune,type = "link")
      lm_prs_tune <- predict(lm_train,data.frame(X_tune),type = "link")
      
      lasso_prs_vad <- predict(lasso_train,X_valid,type = "link")
      ridge_prs_vad <- predict(ridge_train,X_valid,type = "link")
      lm_prs_vad <- predict(lm_train,data.frame(X_valid),type = "link")
      
      
      all_prs_tune <- as.data.frame(cbind(lasso_prs_tune,ridge_prs_tune,lm_prs_tune))
      colnames(all_prs_tune) <- c(paste0("lasso_prs",1:ncol(lasso_prs_tune)),paste0("ridge_prs",1:ncol(ridge_prs_tune)),"lm_prs1")
      all_prs_valid <- as.data.frame(cbind(lasso_prs_vad,ridge_prs_vad,lm_prs_vad))
      colnames(all_prs_valid) <- c(paste0("lasso_prs",1:ncol(lasso_prs_vad)),paste0("ridge_prs",1:ncol(ridge_prs_vad)),"lm_prs1")
      
      mtx <- cor(all_prs_tune)
      drop <- names(all_prs_tune)[apply(mtx,2,function(x){sum(is.na(x))}) == (nrow(mtx) - 1)]
      
      all_prs_tune <- dplyr::select(all_prs_tune, -c(drop))
      all_prs_valid <- dplyr::select(all_prs_valid, -c(drop))
      
      mtx <- cor(all_prs_tune)
      drop <- findCorrelation(mtx,cutoff=0.98)
      drop <- names(all_prs_tune)[drop]
      
      all_prs_tune <- dplyr::select(all_prs_tune, -c(drop))
      all_prs_valid <- dplyr::select(all_prs_valid, -c(drop))
      
      drop <- findLinearCombos(all_prs_tune)$remove
      drop <- names(data.frame(all_prs_tune))[drop]
      
      all_prs_tune <- dplyr::select(all_prs_tune, -c(drop))
      all_prs_valid <- dplyr::select(all_prs_valid, -c(drop))
      
      Results <- Ensemble_Function(x = all_prs_tune,y = pheno_tune[,trait],family = "binary")
      Results$Coefficients[is.na(Results$Coefficients)] <- 0
      
      PRS_Tune <- as.matrix(all_prs_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
      PRS_Validation <- as.matrix(all_prs_valid[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    }
    
    PRS_Tune <- data.frame(IID = pheno_tune$IID,PRS = PRS_Tune)
    PRS_Validation <- data.frame(IID = pheno_validation$IID,PRS = PRS_Validation)
    
    if(trait %in% c("Breast","Prostate")){
      confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }else{
      confounders <- paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
    }
    
    RV_PRS_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),trait)],PRS_Validation)
    RV_PRS_adjusted <- RV_PRS_raw
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    tmp <- data.frame(y = RV_PRS_adjusted[,"PRS"],RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      RV_PRS_adjusted[,"PRS"] <- 0
    }else{
      RV_PRS_adjusted[,"PRS"] <- R/sqrt(y_hat)
    }
    
    RV_PRS_raw_EUR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    RV_PRS_raw_SAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    RV_PRS_raw_AMR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    RV_PRS_raw_AFR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    RV_PRS_raw_EAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    RV_PRS_adjusted_EUR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    RV_PRS_adjusted_SAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    RV_PRS_adjusted_AMR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    RV_PRS_adjusted_AFR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    RV_PRS_adjusted_EAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    RV_PRS_raw_EUR$PRS <- scale(RV_PRS_raw_EUR$PRS)
    RV_PRS_raw_SAS$PRS <- scale(RV_PRS_raw_SAS$PRS)
    RV_PRS_raw_AMR$PRS <- scale(RV_PRS_raw_AMR$PRS)
    RV_PRS_raw_AFR$PRS <- scale(RV_PRS_raw_AFR$PRS)
    RV_PRS_raw_EAS$PRS <- scale(RV_PRS_raw_EAS$PRS)
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2],error = function(e){return(0)})
      return(c(result))
    }
    
    AUC_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
      return(c(result))
    }
    
    beta_validation_raw_EUR <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_EUR,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_raw_EUR, statistic = Beta_Boot, R = 10000)
    beta_raw_EUR_boot <- boot_beta$t
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    
    AUC_validation_raw_EUR <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_EUR[!is.na(RV_PRS_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_raw_EUR, statistic = AUC_Boot, R = 10000)
    AUC_raw_EUR_boot <- boot_AUC$t
    AUC_se_validation_raw_EUR <- sd(boot_AUC$t)
    
    beta_validation_raw_SAS <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_SAS,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_raw_SAS, statistic = Beta_Boot, R = 10000)
    beta_raw_SAS_boot <- boot_beta$t
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    
    AUC_validation_raw_SAS <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_SAS[!is.na(RV_PRS_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_raw_SAS, statistic = AUC_Boot, R = 10000)
    AUC_raw_SAS_boot <- boot_AUC$t
    AUC_se_validation_raw_SAS <- sd(boot_AUC$t)
    
    beta_validation_raw_AMR <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_AMR,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_raw_AMR, statistic = Beta_Boot, R = 10000)
    beta_raw_AMR_boot <- boot_beta$t
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    
    AUC_validation_raw_AMR <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_AMR[!is.na(RV_PRS_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_raw_AMR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AMR_boot <- boot_AUC$t
    AUC_se_validation_raw_AMR <- sd(boot_AUC$t)
    
    beta_validation_raw_AFR <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_AFR,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_raw_AFR, statistic = Beta_Boot, R = 10000)
    beta_raw_AFR_boot <- boot_beta$t
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    
    AUC_validation_raw_AFR <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_AFR[!is.na(RV_PRS_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_raw_AFR, statistic = AUC_Boot, R = 10000)
    AUC_raw_AFR_boot <- boot_AUC$t
    AUC_se_validation_raw_AFR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_EUR <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_EUR,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_adjusted_EUR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EUR_boot <- boot_beta$t
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_EUR <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_EUR[!is.na(RV_PRS_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_adjusted_EUR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_EUR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_EUR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_SAS <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_SAS,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_adjusted_SAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_SAS_boot <- boot_beta$t
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    AUC_validation_adjusted_SAS <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_SAS[!is.na(RV_PRS_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_adjusted_SAS, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_SAS_boot <- boot_AUC$t
    AUC_se_validation_adjusted_SAS <- sd(boot_AUC$t)
    
    beta_validation_adjusted_AMR <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_AMR,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_adjusted_AMR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AMR_boot <- boot_beta$t
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_AMR <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_AMR[!is.na(RV_PRS_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_adjusted_AMR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AMR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_AMR <- sd(boot_AUC$t)
    
    beta_validation_adjusted_AFR <- tryCatch(coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_AFR,family = binomial()))[2],error = function(e){return(0)})
    boot_beta <- boot(data = RV_PRS_adjusted_AFR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AFR_boot <- boot_beta$t
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    AUC_validation_adjusted_AFR <- tryCatch(roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_AFR[!is.na(RV_PRS_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc,error = function(e){return(0)})
    boot_AUC <- boot(data = RV_PRS_adjusted_AFR, statistic = AUC_Boot, R = 10000)
    AUC_adjusted_AFR_boot <- boot_AUC$t
    AUC_se_validation_adjusted_AFR <- sd(boot_AUC$t)
    
    RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR"), 
                                 beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR), 
                                 beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR), 
                                 AUC_raw = c(AUC_validation_raw_EUR,AUC_validation_raw_SAS,AUC_validation_raw_AMR,AUC_validation_raw_AFR),
                                 AUC_se_raw = c(AUC_se_validation_raw_EUR,AUC_se_validation_raw_SAS,AUC_se_validation_raw_AMR,AUC_se_validation_raw_AFR),
                                 beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR), 
                                 beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR), 
                                 AUC_adjusted = c(AUC_validation_adjusted_EUR,AUC_validation_adjusted_SAS,AUC_validation_adjusted_AMR,AUC_validation_adjusted_AFR),
                                 AUC_se_adjusted = c(AUC_se_validation_adjusted_EUR,AUC_se_validation_adjusted_SAS,AUC_se_validation_adjusted_AMR,AUC_se_validation_adjusted_AFR))
    
    RV_Boot_Results <- data.frame(trait = trait,beta_raw_EUR_boot,AUC_raw_EUR_boot,beta_raw_SAS_boot,AUC_raw_SAS_boot,
                                  beta_raw_AMR_boot,AUC_raw_AMR_boot,beta_raw_AFR_boot,AUC_raw_AFR_boot,
                                  beta_adjusted_EUR_boot,AUC_adjusted_EUR_boot,
                                  beta_adjusted_SAS_boot,AUC_adjusted_SAS_boot,beta_adjusted_AMR_boot,AUC_adjusted_AMR_boot,
                                  beta_adjusted_AFR_boot,AUC_adjusted_AFR_boot)
  }else{
    ## Null Models
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_train)
    y_train <- model.null$residual
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
    y_tune <- model.null$residual
    
    model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
    y_validation <- model.null$residual
    
    pheno_train$y_train <- NA
    pheno_train$y_train[!is.na(pheno_train[,trait])] <- y_train
    
    pheno_tune$y_tune <- NA
    pheno_tune$y_tune[!is.na(pheno_tune[,trait])] <- y_tune
    
    pheno_validation$y_validation <- NA
    pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- y_validation
    
    if(ncol(X_train) == 1){
      lm_train <- lm.fit(cbind(1,X_train[!is.na(pheno_train[,trait]),]),pheno_train[!is.na(pheno_train[,trait]),trait])
      lm_train$coefficients[is.na(lm_train$coefficients)] <- 0
      
      lm_prs_tune <- as.numeric(cbind(1,X_tune)%*%matrix(lm_train$coefficients,ncol = 1))
      
      lm_prs_vad <- as.numeric(cbind(1,X_valid)%*%matrix(lm_train$coefficients,ncol = 1))
      
      PRS_Tune <- as.matrix(lm_prs_tune,ncol = 1)
      PRS_Validation <- as.matrix(lm_prs_vad,ncol = 1)
    }else{
      
      lasso_train <- glmnet(X_train[!is.na(pheno_train[,trait]),],pheno_train$y_train[!is.na(pheno_train[,trait])],family = "gaussian",alpha = 1)
      ridge_train <- glmnet(X_train[!is.na(pheno_train[,trait]),],pheno_train$y_train[!is.na(pheno_train[,trait])],family = "gaussian",alpha = 0)
      lm_train <- lm.fit(cbind(1,X_train[!is.na(pheno_train[,trait]),]),pheno_train$y_train[!is.na(pheno_train[,trait])])
      lm_train$coefficients[is.na(lm_train$coefficients)] <- 0
      
      beta_matrix <- as.data.frame(cbind(as.matrix(lasso_train$beta),as.matrix(ridge_train$beta),lm_train$coefficients[-1]))
      colnames(beta_matrix) <- c(paste0("lasso_prs",1:ncol(as.matrix(lasso_train$beta))),paste0("ridge_prs",1:ncol(as.matrix(ridge_train$beta))),paste0("lm_prs",1))
      beta_matrix <- cbind(Train_PVals_All[,c(1:4)],beta_matrix)
      
      
      lasso_prs_tune <- predict(lasso_train,X_tune)
      ridge_prs_tune <- predict(ridge_train,X_tune)
      lm_prs_tune <- as.numeric(cbind(1,X_tune)%*%matrix(lm_train$coefficients,ncol = 1))
      
      lasso_prs_vad <- predict(lasso_train,X_valid)
      ridge_prs_vad <- predict(ridge_train,X_valid)
      lm_prs_vad <- as.numeric(cbind(1,X_valid)%*%matrix(lm_train$coefficients,ncol = 1))
      
      
      all_prs_tune <- as.data.frame(cbind(lasso_prs_tune,ridge_prs_tune,lm_prs_tune))
      colnames(all_prs_tune) <- c(paste0("lasso_prs",1:ncol(lasso_prs_tune)),paste0("ridge_prs",1:ncol(ridge_prs_tune)),"lm_prs1")
      all_prs_valid <- as.data.frame(cbind(lasso_prs_vad,ridge_prs_vad,lm_prs_vad))
      colnames(all_prs_valid) <- c(paste0("lasso_prs",1:ncol(lasso_prs_vad)),paste0("ridge_prs",1:ncol(ridge_prs_vad)),"lm_prs1")
      
      
      mtx <- cor(all_prs_tune)
      drop <- names(all_prs_tune)[apply(mtx,2,function(x){sum(is.na(x))}) == (nrow(mtx) - 1)]
      
      all_prs_tune <- dplyr::select(all_prs_tune, -c(drop))
      all_prs_valid <- dplyr::select(all_prs_valid, -c(drop))
      
      mtx <- cor(all_prs_tune)
      drop <- findCorrelation(mtx,cutoff=0.98)
      drop <- names(all_prs_tune)[drop]
      
      all_prs_tune <- dplyr::select(all_prs_tune, -c(drop))
      all_prs_valid <- dplyr::select(all_prs_valid, -c(drop))
      
      drop <- findLinearCombos(all_prs_tune)$remove
      drop <- names(data.frame(all_prs_tune))[drop]
      
      all_prs_tune <- dplyr::select(all_prs_tune, -c(drop))
      all_prs_valid <- dplyr::select(all_prs_valid, -c(drop))
      
      Results <- Ensemble_Function(x = all_prs_tune,y = pheno_tune[,"y_tune"],family = "continuous")
      Results$Coefficients[is.na(Results$Coefficients)] <- 0
      
      PRS_Tune <- as.matrix(all_prs_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
      PRS_Validation <- as.matrix(all_prs_valid[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
    }
    
    PRS_Tune <- data.frame(IID = pheno_tune$IID,PRS = PRS_Tune)
    PRS_Validation <- data.frame(IID = pheno_validation$IID,PRS = PRS_Validation)
    
    RV_PRS_raw <- inner_join(pheno_validation[,c("IID","age","age2","sex",paste0("pc",1:10),"y_validation")],PRS_Validation)
    RV_PRS_adjusted <- RV_PRS_raw
    
    load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
    
    tmp <- data.frame(y = RV_PRS_adjusted[,"PRS"],RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      RV_PRS_adjusted[,"PRS"] <- 0
    }else{
      RV_PRS_adjusted[,"PRS"] <- R/sqrt(y_hat)
    }
    
    
    RV_PRS_raw_EUR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    RV_PRS_raw_SAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    RV_PRS_raw_AMR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    RV_PRS_raw_AFR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    RV_PRS_raw_EAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    RV_PRS_adjusted_EUR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
    RV_PRS_adjusted_SAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
    RV_PRS_adjusted_AMR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
    RV_PRS_adjusted_AFR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
    RV_PRS_adjusted_EAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
    
    RV_PRS_raw_EUR$y_validation <- scale(RV_PRS_raw_EUR$y_validation)
    RV_PRS_raw_SAS$y_validation <- scale(RV_PRS_raw_SAS$y_validation)
    RV_PRS_raw_AMR$y_validation <- scale(RV_PRS_raw_AMR$y_validation)
    RV_PRS_raw_AFR$y_validation <- scale(RV_PRS_raw_AFR$y_validation)
    RV_PRS_raw_EAS$y_validation <- scale(RV_PRS_raw_EAS$y_validation)
    
    RV_PRS_adjusted_EUR$y_validation <- scale(RV_PRS_adjusted_EUR$y_validation)
    RV_PRS_adjusted_SAS$y_validation <- scale(RV_PRS_adjusted_SAS$y_validation)
    RV_PRS_adjusted_AMR$y_validation <- scale(RV_PRS_adjusted_AMR$y_validation)
    RV_PRS_adjusted_AFR$y_validation <- scale(RV_PRS_adjusted_AFR$y_validation)
    RV_PRS_adjusted_EAS$y_validation <- scale(RV_PRS_adjusted_EAS$y_validation)
    
    RV_PRS_raw_EUR$PRS <- scale(RV_PRS_raw_EUR$PRS)
    RV_PRS_raw_SAS$PRS <- scale(RV_PRS_raw_SAS$PRS)
    RV_PRS_raw_AMR$PRS <- scale(RV_PRS_raw_AMR$PRS)
    RV_PRS_raw_AFR$PRS <- scale(RV_PRS_raw_AFR$PRS)
    RV_PRS_raw_EAS$PRS <- scale(RV_PRS_raw_EAS$PRS)
    
    
    Beta_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = boot_data))[2]
      return(c(result))
    }
    
    R2_Boot <- function(data,indices){
      boot_data <- data[indices, ]
      result <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = boot_data))$r.squared
      return(c(result))
    }
    
    beta_validation_raw_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_EUR))[2]
    boot_beta <- boot(data = RV_PRS_raw_EUR, statistic = Beta_Boot, R = 10000)
    beta_raw_EUR_boot <- boot_beta$t
    beta_se_validation_raw_EUR <- sd(boot_beta$t)
    
    R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_EUR))$r.squared
    boot_R2 <- boot(data = RV_PRS_raw_EUR, statistic = R2_Boot, R = 10000)
    R2_raw_EUR_boot <- boot_R2$t
    R2_se_validation_raw_EUR <- sd(boot_R2$t)
    
    beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_SAS))[2]
    boot_beta <- boot(data = RV_PRS_raw_SAS, statistic = Beta_Boot, R = 10000)
    beta_raw_SAS_boot <- boot_beta$t
    beta_se_validation_raw_SAS <- sd(boot_beta$t)
    
    R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_SAS))$r.squared
    boot_R2 <- boot(data = RV_PRS_raw_SAS, statistic = R2_Boot, R = 10000)
    R2_raw_SAS_boot <- boot_R2$t
    R2_se_validation_raw_SAS <- sd(boot_R2$t)
    
    beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AMR))[2]
    boot_beta <- boot(data = RV_PRS_raw_AMR, statistic = Beta_Boot, R = 10000)
    beta_raw_AMR_boot <- boot_beta$t
    beta_se_validation_raw_AMR <- sd(boot_beta$t)
    
    R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AMR))$r.squared
    boot_R2 <- boot(data = RV_PRS_raw_AMR, statistic = R2_Boot, R = 10000)
    R2_raw_AMR_boot <- boot_R2$t
    R2_se_validation_raw_AMR <- sd(boot_R2$t)
    
    beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AFR))[2]
    boot_beta <- boot(data = RV_PRS_raw_AFR, statistic = Beta_Boot, R = 10000)
    beta_raw_AFR_boot <- boot_beta$t
    beta_se_validation_raw_AFR <- sd(boot_beta$t)
    
    R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AFR))$r.squared
    boot_R2 <- boot(data = RV_PRS_raw_AFR, statistic = R2_Boot, R = 10000)
    R2_raw_AFR_boot <- boot_R2$t
    R2_se_validation_raw_AFR <- sd(boot_R2$t)
    
    beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_EUR))[2]
    boot_beta <- boot(data = RV_PRS_adjusted_EUR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_EUR_boot <- boot_beta$t
    beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
    
    R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_EUR))$r.squared
    boot_R2 <- boot(data = RV_PRS_adjusted_EUR, statistic = R2_Boot, R = 10000)
    R2_adjusted_EUR_boot <- boot_R2$t
    R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
    
    beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_SAS))[2]
    boot_beta <- boot(data = RV_PRS_adjusted_SAS, statistic = Beta_Boot, R = 10000)
    beta_adjusted_SAS_boot <- boot_beta$t
    beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
    
    R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_SAS))$r.squared
    boot_R2 <- boot(data = RV_PRS_adjusted_SAS, statistic = R2_Boot, R = 10000)
    R2_adjusted_SAS_boot <- boot_R2$t
    R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
    
    beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AMR))[2]
    boot_beta <- boot(data = RV_PRS_adjusted_AMR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AMR_boot <- boot_beta$t
    beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AMR))$r.squared
    boot_R2 <- boot(data = RV_PRS_adjusted_AMR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AMR_boot <- boot_R2$t
    R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
    
    beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AFR))[2]
    boot_beta <- boot(data = RV_PRS_adjusted_AFR, statistic = Beta_Boot, R = 10000)
    beta_adjusted_AFR_boot <- boot_beta$t
    beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
    
    R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AFR))$r.squared
    boot_R2 <- boot(data = RV_PRS_adjusted_AFR, statistic = R2_Boot, R = 10000)
    R2_adjusted_AFR_boot <- boot_R2$t
    R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
    
    RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR"), 
                                 beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR), 
                                 beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR), 
                                 R2_raw = c(R2_validation_raw_EUR,R2_validation_raw_SAS,R2_validation_raw_AMR,R2_validation_raw_AFR),
                                 R2_se_raw = c(R2_se_validation_raw_EUR,R2_se_validation_raw_SAS,R2_se_validation_raw_AMR,R2_se_validation_raw_AFR),
                                 beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR), 
                                 beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR), 
                                 R2_adjusted = c(R2_validation_adjusted_EUR,R2_validation_adjusted_SAS,R2_validation_adjusted_AMR,R2_validation_adjusted_AFR),
                                 R2_se_adjusted = c(R2_se_validation_adjusted_EUR,R2_se_validation_adjusted_SAS,R2_se_validation_adjusted_AMR,R2_se_validation_adjusted_AFR))
    
    RV_Boot_Results <- data.frame(trait = trait,beta_raw_EUR_boot,R2_raw_EUR_boot,beta_raw_SAS_boot,R2_raw_SAS_boot,
                                  beta_raw_AMR_boot,R2_raw_AMR_boot,beta_raw_AFR_boot,R2_raw_AFR_boot
                                  ,beta_adjusted_EUR_boot,R2_adjusted_EUR_boot,
                                  beta_adjusted_SAS_boot,R2_adjusted_SAS_boot,beta_adjusted_AMR_boot,R2_adjusted_AMR_boot,
                                  beta_adjusted_AFR_boot,R2_adjusted_AFR_boot)
  } 
}
write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Threshold_Sensitivity_RICE_RV/",trait,"_",thresholds,"Best_Betas.csv"),row.names = FALSE)
write.csv(RV_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Threshold_Sensitivity_RICE_RV/",trait,"_",thresholds,"_Bootstraps.csv"),row.names = FALSE)
