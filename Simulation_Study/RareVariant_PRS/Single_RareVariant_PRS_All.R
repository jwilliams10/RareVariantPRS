rm(list = ls())
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
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

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")

i <- as.numeric(commandArgs(TRUE)[1])

Train_PVals_All <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/coding_sig",i,".Rdata")))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= 1e-03,]
# system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/coding_sig",i,".Rdata")))

## agds dir
agds_dir <- get(load("/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/agds_dir.Rdata"))

## Null Model
obj_nullmodel_train <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Train_Null_Model",i,".RData")))
obj_nullmodel_tune <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Tune_Null_Model",i,".RData")))
obj_nullmodel_validation <- get(load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Validation_Null_Model",i,".RData")))

# system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Train_Null_Model",i,".RData")))
# system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Tune_Null_Model",i,".RData")))
# system(paste0("rm ",paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/nullmodels_staar/Validation_Null_Model",i,".RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22_Annotation_name_catalog.Rdata"))

G_star_gene_centric_coding <- list()

for(j in 1:nrow(Train_PVals_All)){
  ## Chr
  chr <- Train_PVals_All$Chr[j]
  ## Gene name
  gene_name <- Train_PVals_All$Gene[j]
  ## Coding mask
  category <- Train_PVals_All$Category[j]
  
  ### gds file
  gds.path <- agds_dir[chr]
  genofile <- seqOpen(gds.path)
  
  G_star_gene_centric_coding[[j]] <- Gene_Centric_Coding_G_Star(chr=chr,gene_name=gene_name,category=category ,
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
X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_tune$id_include],G_star_gene_centric_coding_tune)
X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_coding_vad)

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")
pheno_train <- Y_train[[i]]
colnames(pheno_train) <- c("IID","Y")
pheno_train <- inner_join(pheno_train,X_train)
X_train <- as.matrix(pheno_train[,3:ncol(pheno_train),drop = FALSE])
model.null <- lm(Y~1,data=pheno_train)
y_train <- model.null$residual

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")
pheno_tune <- Y_tune[[i]]
colnames(pheno_tune) <- c("IID","Y")
pheno_tune <- inner_join(pheno_tune,X_tune)
X_tune <- as.matrix(pheno_tune[,3:ncol(pheno_tune),drop = FALSE])
model.null <- lm(Y~1,data=pheno_tune)
y_tune <- model.null$residual

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")
pheno_valid <- Y_validation[[i]]
colnames(pheno_valid) <- c("IID","Y")
pheno_valid <- inner_join(pheno_valid,X_valid)
X_valid <- as.matrix(pheno_valid[,3:ncol(pheno_valid),drop = FALSE])
model.null <- lm(Y~1,data=pheno_valid)
y_valid <- model.null$residual

if(ncol(X_train) == 1){
  lm_train <- lm.fit(cbind(1,X_train),y_train)
  lm_train$coefficients[is.na(lm_train$coefficients)] <- 0
  
  lm_prs_tune <- as.numeric(cbind(1,X_tune)%*%matrix(lm_train$coefficients,ncol = 1))
  
  lm_prs_vad <- as.numeric(cbind(1,X_valid)%*%matrix(lm_train$coefficients,ncol = 1))
  
  lm_tune_dat <- data.frame(y = y_tune,lm_prs_tune)
  colnames(lm_tune_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_tune_dat) - 1)))
  lm_valid_dat <- data.frame(y = y_valid,lm_prs_vad)
  colnames(lm_valid_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_valid_dat) - 1)))
  
  
  all_prs_tune <- cbind(lm_prs_tune)
  colnames(all_prs_tune) <- c("lm_prs")
  all_prs_valid <- cbind(lm_prs_vad)
  colnames(all_prs_valid) <- c("lm_prs")
  
}else{
  
  lasso_train <- glmnet(X_train,y_train,family = "gaussian",alpha = 1)
  ridge_train <- glmnet(X_train,y_train,family = "gaussian",alpha = 0)
  lm_train <- lm.fit(cbind(1,X_train),y_train)
  lm_train$coefficients[is.na(lm_train$coefficients)] <- 0
  
  lasso_prs_tune <- predict(lasso_train,X_tune)
  ridge_prs_tune <- predict(ridge_train,X_tune)
  lm_prs_tune <- as.numeric(cbind(1,X_tune)%*%matrix(lm_train$coefficients,ncol = 1))
  
  lasso_prs_vad <- predict(lasso_train,X_valid)
  ridge_prs_vad <- predict(ridge_train,X_valid)
  lm_prs_vad <- as.numeric(cbind(1,X_valid)%*%matrix(lm_train$coefficients,ncol = 1))
  
  
  
  
  
  
  lasso_tune_dat <- data.frame(y = y_tune,lasso_prs_tune)
  colnames(lasso_tune_dat) <- c("y",paste0("lasso_prs",1:(ncol(lasso_tune_dat) - 1)))
  lasso_valid_dat <- data.frame(y = y_valid,lasso_prs_vad)
  colnames(lasso_valid_dat) <- c("y",paste0("lasso_prs",1:(ncol(lasso_valid_dat) - 1)))
  
  
  ridge_tune_dat <- data.frame(y = y_tune,ridge_prs_tune)
  colnames(ridge_tune_dat) <- c("y",paste0("ridge_prs",1:(ncol(ridge_tune_dat) - 1)))
  ridge_valid_dat <- data.frame(y = y_valid,ridge_prs_vad)
  colnames(ridge_valid_dat) <- c("y",paste0("ridge_prs",1:(ncol(ridge_valid_dat) - 1)))
  
  lm_tune_dat <- data.frame(y = y_tune,lm_prs_tune)
  colnames(lm_tune_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_tune_dat) - 1)))
  lm_valid_dat <- data.frame(y = y_valid,lm_prs_vad)
  colnames(lm_valid_dat) <- c("y",paste0("lm_prs",1:(ncol(lm_valid_dat) - 1)))
  
  all_prs_tune <- cbind(lasso_prs_tune,ridge_prs_tune,lm_prs_tune)
  colnames(all_prs_tune) <- c(paste0("lasso_prs",1:ncol(lasso_prs_tune)),paste0("ridge_prs",1:ncol(ridge_prs_tune)),"lm_prs")
  all_prs_valid <- cbind(lasso_prs_vad,ridge_prs_vad,lm_prs_vad)
  colnames(all_prs_valid) <- c(paste0("lasso_prs",1:ncol(lasso_prs_vad)),paste0("ridge_prs",1:ncol(ridge_prs_vad)),"lm_prs")
  
}


all_tune <- data.frame(y = y_tune,all_prs_tune)
all_valid <- data.frame(y = y_valid,all_prs_valid)

best_r2 <- vector()
count <- 1
for(j in colnames(all_prs_tune)){
  tmp <- data.frame(y = all_tune$y,x = all_prs_tune[,j])
  best_r2[count] <- summary(lm(y~x,data = tmp))$r.squared
  count <- count + 1
}

best_thresh <- colnames(all_prs_tune)[which.max(best_r2)]
r2_bestoverall_tune <- best_r2[which.max(best_r2)]


all_prs_tune <- as.data.frame(all_prs_tune)
all_prs_valid <- as.data.frame(all_prs_valid)

if(ncol(all_prs_tune)!=1){
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
}

if(ncol(all_prs_tune) == 1){
  tmp_tune <- data.frame(y = y_tune,X = all_prs_tune[,1])
  tmp_valid <- data.frame(X = all_prs_valid[,1])
  mod <- lm(y~X,data = tmp_tune)
  prs_best_tune <- predict(mod,tmp_tune)
  prs_best_vad <- predict(mod,tmp_valid)
}else{
  SL.library <- c(
    "SL.glmnet",
    "SL.ridge",
    "SL.glm",
    "SL.mean"
  )
  
  full_superlearner <- SuperLearner(Y = y_tune, X = all_prs_tune, family = gaussian(),
                                    # For a real analysis we would use V = 10.
                                    # V = 3,
                                    SL.library = SL.library,cvControl = list(V = 20))
  cvsl <- CV.SuperLearner(Y = y_tune, X = all_prs_tune, family = gaussian(),
                          # For a real analysis we would use V = 10.
                          # V = 3,
                          SL.library = SL.library,cvControl = list(V = 20))
  
  best_algorithm <- summary(cvsl)$Table$Algorithm[which.min(summary(cvsl)$Table$Ave)]
  
  a_tune <- predict(full_superlearner, all_prs_tune, onlySL = FALSE)
  a_vad <- predict(full_superlearner, all_prs_valid, onlySL = FALSE)
  
  prs_best_tune_sl <- a_tune$pred
  prs_best_tune_glmnet <- a_tune$library.predict[,1]
  prs_best_tune_ridge <- a_tune$library.predict[,2]
  prs_best_tune_glm <- a_tune$library.predict[,3]
  prs_best_tune_mean <- a_tune$library.predict[,4]
  
  prs_best_vad_sl <- a_vad$pred
  prs_best_vad_glmnet <- a_vad$library.predict[,1]
  prs_best_vad_ridge <- a_vad$library.predict[,2]
  prs_best_vad_glm <- a_vad$library.predict[,3]
  prs_best_vad_mean <- a_vad$library.predict[,4]
  
  if(best_algorithm == "SL.glmnet_All"){
    #final
    prs_best_tune <- prs_best_tune_glmnet
    prs_best_vad <- prs_best_vad_glmnet
  }else if(best_algorithm == "SL.ridge_All"){
    #final
    prs_best_tune <- prs_best_tune_ridge
    prs_best_vad <- prs_best_vad_ridge
  }else if(best_algorithm == "SL.glm_All"){
    #final
    prs_best_tune <- prs_best_tune_glm
    prs_best_vad <- prs_best_vad_glm
  }else if(best_algorithm == "SL.mean_All"){
    #final
    prs_best_tune <- prs_best_tune_mean
    prs_best_vad <- prs_best_vad_mean
  } else {
    prs_best_tune <- prs_best_tune_sl
    prs_best_vad <- prs_best_vad_sl
  }
}

tune_dat_sl_R2 <- data.frame(y = y_tune,x = prs_best_tune)
valid_dat_sl_R2 <- data.frame(y = y_valid,x = prs_best_vad)

r2_sl_tune <- summary(lm(y~x,data = tune_dat_sl_R2))$r.squared

if(r2_sl_tune < r2_bestoverall_tune){
  prs_best_tune_sl <- all_tune[,best_thresh]
  prs_best_vad_sl <- all_valid[,best_thresh] 
  
  tune_dat_sl_R2 <- data.frame(y = y_tune,x = prs_best_tune_sl)
  valid_dat_sl_R2 <- data.frame(y = y_valid,x = prs_best_vad_sl)
}





load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

RV_PRS <- data.frame(IID = pheno_valid$IID,Y = valid_dat_sl_R2[,1],RV_PRS = valid_dat_sl_R2[,2])

write.csv(RV_PRS,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/BestPRS",i,".csv"),row.names = FALSE)

RV_PRS_raw <- RV_PRS
RV_PRS_adjusted <- RV_PRS
RV_PRS_adjusted <- inner_join(RV_PRS_adjusted,ukb_pheno[,c("IID","pc1","pc2","pc3","pc4","pc5")])


tmp <- data.frame(y = RV_PRS_adjusted[,"RV_PRS"],RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
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
  RV_PRS_adjusted[,"RV_PRS"] <- 0
}else{
  RV_PRS_adjusted[,"RV_PRS"] <- R/sqrt(y_hat)
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

RV_PRS_raw_EUR$Y <- scale(RV_PRS_raw_EUR$Y)
RV_PRS_raw_SAS$Y <- scale(RV_PRS_raw_SAS$Y)
RV_PRS_raw_AMR$Y <- scale(RV_PRS_raw_AMR$Y)
RV_PRS_raw_AFR$Y <- scale(RV_PRS_raw_AFR$Y)
RV_PRS_raw_EAS$Y <- scale(RV_PRS_raw_EAS$Y)

RV_PRS_adjusted_EUR$Y <- scale(RV_PRS_adjusted_EUR$Y)
RV_PRS_adjusted_SAS$Y <- scale(RV_PRS_adjusted_SAS$Y)
RV_PRS_adjusted_AMR$Y <- scale(RV_PRS_adjusted_AMR$Y)
RV_PRS_adjusted_AFR$Y <- scale(RV_PRS_adjusted_AFR$Y)
RV_PRS_adjusted_EAS$Y <- scale(RV_PRS_adjusted_EAS$Y)

RV_PRS_raw_EUR$RV_PRS <- scale(RV_PRS_raw_EUR$RV_PRS)
RV_PRS_raw_SAS$RV_PRS <- scale(RV_PRS_raw_SAS$RV_PRS)
RV_PRS_raw_AMR$RV_PRS <- scale(RV_PRS_raw_AMR$RV_PRS)
RV_PRS_raw_AFR$RV_PRS <- scale(RV_PRS_raw_AFR$RV_PRS)
RV_PRS_raw_EAS$RV_PRS <- scale(RV_PRS_raw_EAS$RV_PRS)


best_beta_raw_RV_EUR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_EUR))[2]
se_beta_raw_RV_EUR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_EUR))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_SAS))[2]
se_beta_raw_RV_SAS <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_SAS))$coefficients[2,2]
best_beta_raw_RV_AMR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_AMR))[2]
se_beta_raw_RV_AMR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_AMR))$coefficients[2,2]
best_beta_raw_RV_AFR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_AFR))[2]
se_beta_raw_RV_AFR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_AFR))$coefficients[2,2]
best_beta_raw_RV_EAS <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_EAS))[2]
se_beta_raw_RV_EAS <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_EAS))$coefficients[2,2]

best_beta_adjusted_RV_EUR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_EUR))[2]
se_beta_adjusted_RV_EUR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_EUR))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_SAS))[2]
se_beta_adjusted_RV_SAS <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_SAS))$coefficients[2,2]
best_beta_adjusted_RV_AMR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_AMR))[2]
se_beta_adjusted_RV_AMR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_AMR))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_AFR))[2]
se_beta_adjusted_RV_AFR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_AFR))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_EAS))[2]
se_beta_adjusted_RV_EAS <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_EAS))$coefficients[2,2]

RV_PRS_Results <- data.frame(i = i,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_SAS,best_beta_raw_RV_AMR,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_SAS,se_beta_raw_RV_AMR,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_AMR,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                             se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_AMR,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS))

write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/Best_Betas",i,".csv"),row.names = FALSE)
