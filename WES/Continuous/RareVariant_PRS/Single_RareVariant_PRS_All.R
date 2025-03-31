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

source("~/RareVariantPRS/Simulation_Study/FullData/G_star/g_star_gene_centric_coding.R")

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "BMI"
}else if(trait == 2){
  trait <- "TC"
}else if(trait == 3){
  trait <- "HDL"
}else if(trait == 4){
  trait <- "LDL"
}else if(trait == 5){
  trait <- "logTG"
}else{
  trait <- "Height"
}

Train_PVals_All <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= 1e-03,]

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))

## Null Model
obj_nullmodel_train <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Train_Null_Model.RData")))
obj_nullmodel_tune <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/nullmodels_staar/",trait,"_Validation_Null_Model.RData")))

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
rm(G_star_gene_centric_coding)


X_train <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_train$id_include],G_star_gene_centric_coding_train)
pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
pheno_train <- inner_join(pheno_train,X_train)
X_train <- as.matrix(pheno_train[,27:ncol(pheno_train),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_train)
pheno_train$y_train <- NA
pheno_train$y_train[!is.na(pheno_train[,trait])] <- model.null$residual

X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_tune$id_include],G_star_gene_centric_coding_tune)
pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tune <- inner_join(pheno_tune,X_tune)
X_tune <- as.matrix(pheno_tune[,27:ncol(pheno_tune),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
pheno_tune$y_tune <- NA
pheno_tune$y_tune[!is.na(pheno_tune[,trait])] <- model.null$residual

X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_coding_vad)
pheno_valid <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_valid <- inner_join(pheno_valid,X_valid)
X_valid <- as.matrix(pheno_valid[,27:ncol(pheno_valid),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_valid)
pheno_valid$y_validation <- NA
pheno_valid$y_validation[!is.na(pheno_valid[,trait])] <- model.null$residual


lasso_train <- glmnet(X_train,pheno_train$y_train,family = "gaussian",alpha = 1)
ridge_train <- glmnet(X_train,pheno_train$y_train,family = "gaussian",alpha = 0)
lm_train <- lm.fit(cbind(1,X_train),pheno_train$y_train)
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


lasso_tune_dat <- data.frame(y = pheno_tune$y_tune,lasso_prs_tune)
colnames(lasso_tune_dat) <- c("y",paste0("lasso_prs",1:(ncol(lasso_tune_dat) - 1)))
lasso_valid_dat <- data.frame(y = pheno_valid$y_validation,lasso_prs_vad)
colnames(lasso_valid_dat) <- c("y",paste0("lasso_prs",1:(ncol(lasso_valid_dat) - 1)))

ridge_tune_dat <- data.frame(y = pheno_tune$y_tune,ridge_prs_tune)
colnames(ridge_tune_dat) <- c("y",paste0("ridge_prs",1:(ncol(ridge_tune_dat) - 1)))
ridge_valid_dat <- data.frame(y = pheno_valid$y_validation ,ridge_prs_vad)
colnames(ridge_valid_dat) <- c("y",paste0("ridge_prs",1:(ncol(ridge_valid_dat) - 1)))

lm_tune_dat <- data.frame(y = pheno_tune$y_tune,lm_prs_tune)
colnames(lm_tune_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_tune_dat) - 1)))
lm_valid_dat <- data.frame(y = pheno_valid$y_validation,lm_prs_vad)
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
  
  lasso_mod <- cv.glmnet(x,y,family = "gaussian",alpha = 1,type.measure = "mse",nfold = 10)
  ridge_mod <- cv.glmnet(x,y,family = "gaussian",alpha = 0,type.measure = "mse",nfold = 10)
  
  lasso_prediction_x <- predict(lasso_mod, x)
  ridge_prediction_x <- predict(ridge_mod, x)
  
  ensemble_mod <- lm(y~.,data = data.frame(lasso_prediction_x,ridge_prediction_x))
  
  ensemble_prediction_x <- ensemble_mod$fitted
  
  coefficients_x <- coef(lm(y~.,data.frame(y = ensemble_prediction_x,x)))
  return(list(Coefficients = coefficients_x))
}
Ensemble_Function_Binary<- function(x,y){
  x <- as.matrix(x[!is.na(y),])
  y <- y[!is.na(y)]
  
  lasso_mod <- cv.glmnet(x,y,family = "binomial",alpha = 1,type.measure = "auc")
  ridge_mod <- cv.glmnet(x,y,family = "binomial",alpha = 0,type.measure = "auc")
  
  lasso_prediction_x <- predict(lasso_mod, x,type = "link")
  ridge_prediction_x <- predict(ridge_mod, x,type = "link")
  
  ensemble_mod <- glm(y~.,data = data.frame(lasso_prediction_x,ridge_prediction_x),family = binomial())
  ensemble_prediction_x <- predict(ensemble_mod,data.frame(lasso_prediction_x,ridge_prediction_x),type = "link")
  
  coefficients_x <- coef(lm(y~.,data.frame(y = ensemble_prediction_x,x)))
  return(list(Coefficients = coefficients_x))
}
Ensemble_Function <- function(x,y,family = c("continuous","binary")){
  if(family == "continuous"){
    return(Ensemble_Function_Continuous(x,y))
  }else{
    return(Ensemble_Function_Binary(x,y))
  }
}

Results <- Ensemble_Function(x = all_prs_tune,y = pheno_tune[,"y_tune"],family = "continuous")
Results$Coefficients[is.na(Results$Coefficients)] <- 0
Final_Coefficients <- data.frame(beta_matrix[,1:4],BETA = as.matrix(beta_matrix[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1))
write.csv(Final_Coefficients,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"),row.names = FALSE)
PRS_Tune <- as.matrix(all_prs_tune[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)
PRS_Validation <- as.matrix(all_prs_valid[,names(Results$Coefficients)[-1]]) %*% matrix(Results$Coefficients[-1],ncol = 1)

PRS_Tune <- data.frame(IID = pheno_tune$IID, PRS = PRS_Tune)
PRS_Validation <- data.frame(IID = pheno_valid$IID, PRS = PRS_Validation)

write.csv(PRS_Tune,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_BestPRS_Tune.csv"),row.names = FALSE)
write.csv(PRS_Validation,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_BestPRS_Validation.csv"),row.names = FALSE)

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

RV_PRS <- inner_join(pheno_valid[,c("IID","y_validation","age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")],PRS_Validation)

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
boot_beta <- boot(data = RV_PRS_raw_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EUR <- sd(boot_beta$t)
beta_lower_validation_raw_EUR <- beta_ci$basic[4]
beta_upper_validation_raw_EUR <- beta_ci$basic[5]

R2_validation_raw_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_EUR))$r.squared
boot_R2 <- boot(data = RV_PRS_raw_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EUR <- sd(boot_R2$t)
R2_lower_validation_raw_EUR <- R2_ci$basic[4]
R2_upper_validation_raw_EUR <- R2_ci$basic[5]

beta_validation_raw_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_SAS))[2]
boot_beta <- boot(data = RV_PRS_raw_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_SAS <- sd(boot_beta$t)
beta_lower_validation_raw_SAS <- beta_ci$basic[4]
beta_upper_validation_raw_SAS <- beta_ci$basic[5]

R2_validation_raw_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_SAS))$r.squared
boot_R2 <- boot(data = RV_PRS_raw_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_SAS <- sd(boot_R2$t)
R2_lower_validation_raw_SAS <- R2_ci$basic[4]
R2_upper_validation_raw_SAS <- R2_ci$basic[5]

beta_validation_raw_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AMR))[2]
boot_beta <- boot(data = RV_PRS_raw_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AMR <- sd(boot_beta$t)
beta_lower_validation_raw_AMR <- beta_ci$basic[4]
beta_upper_validation_raw_AMR <- beta_ci$basic[5]

R2_validation_raw_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AMR))$r.squared
boot_R2 <- boot(data = RV_PRS_raw_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AMR <- sd(boot_R2$t)
R2_lower_validation_raw_AMR <- R2_ci$basic[4]
R2_upper_validation_raw_AMR <- R2_ci$basic[5]

beta_validation_raw_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AFR))[2]
boot_beta <- boot(data = RV_PRS_raw_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_AFR <- sd(boot_beta$t)
beta_lower_validation_raw_AFR <- beta_ci$basic[4]
beta_upper_validation_raw_AFR <- beta_ci$basic[5]

R2_validation_raw_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_AFR))$r.squared
boot_R2 <- boot(data = RV_PRS_raw_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_AFR <- sd(boot_R2$t)
R2_lower_validation_raw_AFR <- R2_ci$basic[4]
R2_upper_validation_raw_AFR <- R2_ci$basic[5]

beta_validation_raw_EAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_EAS))[2]
boot_beta <- boot(data = RV_PRS_raw_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_raw_EAS <- sd(boot_beta$t)
beta_lower_validation_raw_EAS <- beta_ci$basic[4]
beta_upper_validation_raw_EAS <- beta_ci$basic[5]

R2_validation_raw_EAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_raw_EAS))$r.squared
boot_R2 <- boot(data = RV_PRS_raw_EAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_raw_EAS <- sd(boot_R2$t)
R2_lower_validation_raw_EAS <- R2_ci$basic[4]
R2_upper_validation_raw_EAS <- R2_ci$basic[5]

beta_validation_adjusted_EUR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_EUR))[2]
boot_beta <- boot(data = RV_PRS_adjusted_EUR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EUR <- sd(boot_beta$t)
beta_lower_validation_adjusted_EUR <- beta_ci$basic[4]
beta_upper_validation_adjusted_EUR <- beta_ci$basic[5]

R2_validation_adjusted_EUR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_EUR))$r.squared
boot_R2 <- boot(data = RV_PRS_adjusted_EUR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EUR <- sd(boot_R2$t)
R2_lower_validation_adjusted_EUR <- R2_ci$basic[4]
R2_upper_validation_adjusted_EUR <- R2_ci$basic[5]

beta_validation_adjusted_SAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_SAS))[2]
boot_beta <- boot(data = RV_PRS_adjusted_SAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_SAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_SAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_SAS <- beta_ci$basic[5]

R2_validation_adjusted_SAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_SAS))$r.squared
boot_R2 <- boot(data = RV_PRS_adjusted_SAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_SAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_SAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_SAS <- R2_ci$basic[5]

beta_validation_adjusted_AMR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AMR))[2]
boot_beta <- boot(data = RV_PRS_adjusted_AMR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AMR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AMR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AMR <- beta_ci$basic[5]

R2_validation_adjusted_AMR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AMR))$r.squared
boot_R2 <- boot(data = RV_PRS_adjusted_AMR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AMR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AMR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AMR <- R2_ci$basic[5]

beta_validation_adjusted_AFR <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AFR))[2]
boot_beta <- boot(data = RV_PRS_adjusted_AFR, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_AFR <- sd(boot_beta$t)
beta_lower_validation_adjusted_AFR <- beta_ci$basic[4]
beta_upper_validation_adjusted_AFR <- beta_ci$basic[5]

R2_validation_adjusted_AFR <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_AFR))$r.squared
boot_R2 <- boot(data = RV_PRS_adjusted_AFR, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_AFR <- sd(boot_R2$t)
R2_lower_validation_adjusted_AFR <- R2_ci$basic[4]
R2_upper_validation_adjusted_AFR <- R2_ci$basic[5]

beta_validation_adjusted_EAS <- coef(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_EAS))[2]
boot_beta <- boot(data = RV_PRS_adjusted_EAS, statistic = Beta_Boot, R = 1000)
beta_ci <- boot.ci(boot_beta, type = "basic")
beta_se_validation_adjusted_EAS <- sd(boot_beta$t)
beta_lower_validation_adjusted_EAS <- beta_ci$basic[4]
beta_upper_validation_adjusted_EAS <- beta_ci$basic[5]

R2_validation_adjusted_EAS <- summary(lm(as.formula(paste0("y_validation~","PRS")),data = RV_PRS_adjusted_EAS))$r.squared
boot_R2 <- boot(data = RV_PRS_adjusted_EAS, statistic = R2_Boot, R = 1000)
R2_ci <- boot.ci(boot_R2, type = "basic")
R2_se_validation_adjusted_EAS <- sd(boot_R2$t)
R2_lower_validation_adjusted_EAS <- R2_ci$basic[4]
R2_upper_validation_adjusted_EAS <- R2_ci$basic[5]

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
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

write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"Best_Betas.csv"),row.names = FALSE)
