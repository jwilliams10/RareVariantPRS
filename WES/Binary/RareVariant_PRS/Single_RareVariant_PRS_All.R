rm(list = ls())
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS,lib.loc = "/usr/local/apps/R/4.3/site-library_4.3.2")
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readr)
library(dplyr)
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(stringr)
library(glmnet)
library(RISCA)

source("~/RareVariantPRS/Simulation_Study/FullData/G_star/g_star_gene_centric_coding.R")

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "Asthma"
}else if(trait == 2){
  trait <- "CAD"
}else if(trait == 3){
  trait <- "T2D"
}else if(trait == 4){
  trait <- "Breast"
}else{
  trait <- "Prostate"
}

if(trait %in% c("Breast","Prostate")){
  confounders <- paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}else{
  confounders <-  paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")
}

Train_PVals_All <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= 1e-03,]

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))

## Null Model
obj_nullmodel_train <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Train_Null_Model.RData")))
obj_nullmodel_tune <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Validation_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)

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

col_remove <- apply(G_star_gene_centric_coding,2,function(x){sum(x != 0)}) > 10 & colSums(G_star_gene_centric_coding) > 10 
G_star_gene_centric_coding <- G_star_gene_centric_coding[,col_remove,drop = FALSE]

Train_PVals_All <- Train_PVals_All[col_remove,]

ids_gstar <- obj_nullmodel$id_include

G_star_gene_centric_coding_train <- G_star_gene_centric_coding[ids_gstar %in% obj_nullmodel_train$id_include,]
G_star_gene_centric_coding_tune <- G_star_gene_centric_coding[ids_gstar %in% obj_nullmodel_tune$id_include,]
G_star_gene_centric_coding_vad <- G_star_gene_centric_coding[ids_gstar %in% obj_nullmodel_validation$id_include,]


X_train <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_train$id_include],G_star_gene_centric_coding_train)
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- inner_join(pheno_train,X_train)
X_train <- as.matrix(pheno_train[,27:ncol(pheno_train),drop = FALSE])

X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_tune$id_include],G_star_gene_centric_coding_tune)
pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tune <- inner_join(pheno_tune,X_tune)
X_tune <- as.matrix(pheno_tune[,27:ncol(pheno_tune),drop = FALSE])

X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_coding_vad)
pheno_valid <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_valid <- inner_join(pheno_valid,X_valid)
X_valid <- as.matrix(pheno_valid[,27:ncol(pheno_valid),drop = FALSE])


lasso_train <- glmnet(X_train,pheno_train[,trait],family = "binomial",alpha = 1)
ridge_train <- glmnet(X_train,pheno_train[,trait],family = "binomial",alpha = 0)
lm_train <- glm(y~.,data = data.frame(y = pheno_train[,trait], X_train),family = binomial())

beta_matrix <- as.data.frame(cbind(as.matrix(lasso_train$beta),as.matrix(ridge_train$beta),lm_train$coefficients[-1]))
colnames(beta_matrix) <- c(paste0("lasso_prs",1:ncol(as.matrix(lasso_train$beta))),paste0("ridge_prs",1:ncol(as.matrix(ridge_train$beta))),paste0("lm_prs",1))
beta_matrix <- cbind(Train_PVals_All[,c(1:4)],beta_matrix)


lasso_prs_tune <- predict(lasso_train,X_tune,type = "link")
ridge_prs_tune <- predict(ridge_train,X_tune,type = "link")
lm_prs_tune <- predict(lm_train,data.frame(X_tune),type = "link")

lasso_prs_vad <- predict(lasso_train,X_valid,type = "link")
ridge_prs_vad <- predict(ridge_train,X_valid,type = "link")
lm_prs_vad <- predict(lm_train,data.frame(X_valid),type = "link")

lasso_tune_dat <- data.frame(y = pheno_tune[,trait],lasso_prs_tune)
colnames(lasso_tune_dat) <- c("y",paste0("lasso_prs",1:(ncol(lasso_tune_dat) - 1)))
lasso_valid_dat <- data.frame(y = pheno_valid[,trait],lasso_prs_vad)
colnames(lasso_valid_dat) <- c("y",paste0("lasso_prs",1:(ncol(lasso_valid_dat) - 1)))

ridge_tune_dat <- data.frame(y = pheno_tune[,trait],ridge_prs_tune)
colnames(ridge_tune_dat) <- c("y",paste0("ridge_prs",1:(ncol(ridge_tune_dat) - 1)))
ridge_valid_dat <- data.frame(y = pheno_valid[,trait],ridge_prs_vad)
colnames(ridge_valid_dat) <- c("y",paste0("ridge_prs",1:(ncol(ridge_valid_dat) - 1)))

lm_tune_dat <- data.frame(y = pheno_tune[,trait],lm_prs_tune)
colnames(lm_tune_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_tune_dat) - 1)))
lm_valid_dat <- data.frame(y = pheno_valid[,trait],lm_prs_vad)
colnames(lm_valid_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_valid_dat) - 1)))


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

Results <- Ensemble_Function(x = all_prs_tune,y = pheno_tune[,trait],family = "binary")
Results$Coefficients[is.na(Results$Coefficients)] <- 0
Final_Coefficients <- data.frame(beta_matrix[,1:4],BETA = as.matrix(beta_matrix[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1))
write.csv(Final_Coefficients,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"),row.names = FALSE)
PRS_Tune <- as.matrix(all_prs_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Validation <- as.matrix(all_prs_valid[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)

PRS_Tune <- data.frame(IID = pheno_tune$IID, PRS = PRS_Tune)
PRS_Validation <- data.frame(IID = pheno_valid$IID, PRS = PRS_Validation)

write.csv(PRS_Tune,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_BestPRS_Tune.csv"),row.names = FALSE)
write.csv(PRS_Validation,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_BestPRS_Validation.csv"),row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

RV_PRS <- inner_join(pheno_valid[,c("IID",trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],PRS_Validation)

RV_PRS_raw <- RV_PRS
RV_PRS_adjusted <- RV_PRS

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
  result <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = boot_data,family = binomial()))[2]
  return(c(result))
}

AUC_Boot <- function(data,indices){
  boot_data <- data[indices, ]
  result <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = boot_data[!is.na(boot_data[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
  return(c(result))
}

beta_validation_raw_EUR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_EUR,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_raw_EUR, statistic = Beta_Boot, R = 10000)
beta_raw_EUR_boot <- boot_beta$t
beta_se_validation_raw_EUR <- sd(boot_beta$t)

AUC_validation_raw_EUR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_EUR[!is.na(RV_PRS_raw_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_raw_EUR, statistic = AUC_Boot, R = 10000)
AUC_raw_EUR_boot <- boot_AUC$t
AUC_se_validation_raw_EUR <- sd(boot_AUC$t)

beta_validation_raw_SAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_SAS,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_raw_SAS, statistic = Beta_Boot, R = 10000)
beta_raw_SAS_boot <- boot_beta$t
beta_se_validation_raw_SAS <- sd(boot_beta$t)

AUC_validation_raw_SAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_SAS[!is.na(RV_PRS_raw_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_raw_SAS, statistic = AUC_Boot, R = 10000)
AUC_raw_SAS_boot <- boot_AUC$t
AUC_se_validation_raw_SAS <- sd(boot_AUC$t)

beta_validation_raw_AMR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_AMR,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_raw_AMR, statistic = Beta_Boot, R = 10000)
beta_raw_AMR_boot <- boot_beta$t
beta_se_validation_raw_AMR <- sd(boot_beta$t)

AUC_validation_raw_AMR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_AMR[!is.na(RV_PRS_raw_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_raw_AMR, statistic = AUC_Boot, R = 10000)
AUC_raw_AMR_boot <- boot_AUC$t
AUC_se_validation_raw_AMR <- sd(boot_AUC$t)

beta_validation_raw_AFR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_AFR,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_raw_AFR, statistic = Beta_Boot, R = 10000)
beta_raw_AFR_boot <- boot_beta$t
beta_se_validation_raw_AFR <- sd(boot_beta$t)

AUC_validation_raw_AFR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_AFR[!is.na(RV_PRS_raw_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_raw_AFR, statistic = AUC_Boot, R = 10000)
AUC_raw_AFR_boot <- boot_AUC$t
AUC_se_validation_raw_AFR <- sd(boot_AUC$t)

beta_validation_raw_EAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_raw_EAS,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_raw_EAS, statistic = Beta_Boot, R = 10000)
beta_raw_EAS_boot <- boot_beta$t
beta_se_validation_raw_EAS <- sd(boot_beta$t)

AUC_validation_raw_EAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_raw_EAS[!is.na(RV_PRS_raw_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_raw_EAS, statistic = AUC_Boot, R = 10000)
AUC_raw_EAS_boot <- boot_AUC$t
AUC_se_validation_raw_EAS <- sd(boot_AUC$t)

beta_validation_adjusted_EUR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_EUR,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_adjusted_EUR, statistic = Beta_Boot, R = 10000)
beta_adjusted_EUR_boot <- boot_beta$t
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)

AUC_validation_adjusted_EUR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_EUR[!is.na(RV_PRS_adjusted_EUR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_adjusted_EUR, statistic = AUC_Boot, R = 10000)
AUC_adjusted_EUR_boot <- boot_AUC$t
AUC_se_validation_adjusted_EUR <- sd(boot_AUC$t)

beta_validation_adjusted_SAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_SAS,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_adjusted_SAS, statistic = Beta_Boot, R = 10000)
beta_adjusted_SAS_boot <- boot_beta$t
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)

AUC_validation_adjusted_SAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_SAS[!is.na(RV_PRS_adjusted_SAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_adjusted_SAS, statistic = AUC_Boot, R = 10000)
AUC_adjusted_SAS_boot <- boot_AUC$t
AUC_se_validation_adjusted_SAS <- sd(boot_AUC$t)

beta_validation_adjusted_AMR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_AMR,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_adjusted_AMR, statistic = Beta_Boot, R = 10000)
beta_adjusted_AMR_boot <- boot_beta$t
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)

AUC_validation_adjusted_AMR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_AMR[!is.na(RV_PRS_adjusted_AMR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_adjusted_AMR, statistic = AUC_Boot, R = 10000)
AUC_adjusted_AMR_boot <- boot_AUC$t
AUC_se_validation_adjusted_AMR <- sd(boot_AUC$t)

beta_validation_adjusted_AFR <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_AFR,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_adjusted_AFR, statistic = Beta_Boot, R = 10000)
beta_adjusted_AFR_boot <- boot_beta$t
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)

AUC_validation_adjusted_AFR <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_AFR[!is.na(RV_PRS_adjusted_AFR[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_adjusted_AFR, statistic = AUC_Boot, R = 10000)
AUC_adjusted_AFR_boot <- boot_AUC$t
AUC_se_validation_adjusted_AFR <- sd(boot_AUC$t)

beta_validation_adjusted_EAS <- coef(glm(as.formula(paste0(trait,"~","PRS","+",gsub("~","",confounders))),data = RV_PRS_adjusted_EAS,family = binomial()))[2]
boot_beta <- boot(data = RV_PRS_adjusted_EAS, statistic = Beta_Boot, R = 10000)
beta_adjusted_EAS_boot <- boot_beta$t
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)

AUC_validation_adjusted_EAS <- roc.binary(status = trait,variable = "PRS",confounders = as.formula(confounders),data = RV_PRS_adjusted_EAS[!is.na(RV_PRS_adjusted_EAS[,trait]),],precision=seq(0.05,0.95, by=0.05))$auc
boot_AUC <- boot(data = RV_PRS_adjusted_EAS, statistic = AUC_Boot, R = 10000)
AUC_adjusted_EAS_boot <- boot_AUC$t
AUC_se_validation_adjusted_EAS <- sd(boot_AUC$t)

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(beta_validation_raw_EUR,beta_validation_raw_SAS,beta_validation_raw_AMR,beta_validation_raw_AFR,beta_validation_raw_EAS), 
                             beta_se_raw = c(beta_se_validation_raw_EUR,beta_se_validation_raw_SAS,beta_se_validation_raw_AMR,beta_se_validation_raw_AFR,beta_se_validation_raw_EAS), 
                             AUC_raw = c(AUC_validation_raw_EUR,AUC_validation_raw_SAS,AUC_validation_raw_AMR,AUC_validation_raw_AFR,AUC_validation_raw_EAS),
                             AUC_se_raw = c(AUC_se_validation_raw_EUR,AUC_se_validation_raw_SAS,AUC_se_validation_raw_AMR,AUC_se_validation_raw_AFR,AUC_se_validation_raw_EAS),
                             beta_adjusted = c(beta_validation_adjusted_EUR,beta_validation_adjusted_SAS,beta_validation_adjusted_AMR,beta_validation_adjusted_AFR,beta_validation_adjusted_EAS), 
                             beta_se_adjusted = c(beta_se_validation_adjusted_EUR,beta_se_validation_adjusted_SAS,beta_se_validation_adjusted_AMR,beta_se_validation_adjusted_AFR,beta_se_validation_adjusted_EAS), 
                             AUC_adjusted = c(AUC_validation_adjusted_EUR,AUC_validation_adjusted_SAS,AUC_validation_adjusted_AMR,AUC_validation_adjusted_AFR,AUC_validation_adjusted_EAS),
                             AUC_se_adjusted = c(AUC_se_validation_adjusted_EUR,AUC_se_validation_adjusted_SAS,AUC_se_validation_adjusted_AMR,AUC_se_validation_adjusted_AFR,AUC_se_validation_adjusted_EAS))

RV_Boot_Results <- data.frame(trait = trait,beta_raw_EUR_boot,AUC_raw_EUR_boot,beta_raw_SAS_boot,AUC_raw_SAS_boot,
                              beta_raw_AMR_boot,AUC_raw_AMR_boot,beta_raw_AFR_boot,AUC_raw_AFR_boot,
                              beta_raw_EAS_boot,AUC_raw_EAS_boot,beta_adjusted_EUR_boot,AUC_adjusted_EUR_boot,
                              beta_adjusted_SAS_boot,AUC_adjusted_SAS_boot,beta_adjusted_AMR_boot,AUC_adjusted_AMR_boot,
                              beta_adjusted_AFR_boot,AUC_adjusted_AFR_boot,beta_adjusted_EAS_boot,AUC_adjusted_EAS_boot)

write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"Best_Betas.csv"),row.names = FALSE)
write.csv(RV_Boot_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_RareVariants_PRS/",trait,"_Bootstraps.csv"),row.names = FALSE) 