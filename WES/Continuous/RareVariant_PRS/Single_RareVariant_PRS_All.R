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
y_train <- model.null$residual

X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_tune$id_include],G_star_gene_centric_coding_tune)
pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
pheno_tune <- inner_join(pheno_tune,X_tune)
X_tune <- as.matrix(pheno_tune[,27:ncol(pheno_tune),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
y_tune <- model.null$residual

X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_coding_vad)
pheno_valid <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
pheno_valid <- inner_join(pheno_valid,X_valid)
X_valid <- as.matrix(pheno_valid[,27:ncol(pheno_valid),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_valid)
y_valid <- model.null$residual

lasso_train <- glmnet(X_train,y_train,family = "gaussian",alpha = 1)
ridge_train <- glmnet(X_train,y_train,family = "gaussian",alpha = 0)
lm_train <- lm.fit(cbind(1,X_train),y_train)
lm_train$coefficients[is.na(lm_train$coefficients)] <- 0

beta_matrix <- cbind(as.matrix(lasso_train$beta),as.matrix(ridge_train$beta),lm_train$coefficients[-1])
beta_matrix <- as.data.frame(beta_matrix)
colnames(beta_matrix) <- c(paste0("lasso_prs",1:ncol(as.matrix(lasso_train$beta))),paste0("ridge_prs",1:ncol(as.matrix(ridge_train$beta))),paste0("lm_prs",1))

beta_matrix <- cbind(Train_PVals_All[,c(1:4)],beta_matrix)


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
colnames(all_prs_tune) <- c(paste0("lasso_prs",1:ncol(lasso_prs_tune)),paste0("ridge_prs",1:ncol(ridge_prs_tune)),"lm_prs1")
all_prs_valid <- cbind(lasso_prs_vad,ridge_prs_vad,lm_prs_vad)
colnames(all_prs_valid) <- c(paste0("lasso_prs",1:ncol(lasso_prs_vad)),paste0("ridge_prs",1:ncol(ridge_prs_vad)),"lm_prs1")

all_tune <- data.frame(y = y_tune,all_prs_tune)
all_valid <- data.frame(y = y_valid,all_prs_valid)

best_r2 <- vector()
count <- 1
for(i in colnames(all_prs_tune)){
  tmp <- data.frame(y = all_tune$y,x = all_prs_tune[,i])
  best_r2[count] <- summary(lm(y~x,data = tmp))$r.squared
  count <- count + 1
}

best_thresh <- colnames(all_prs_tune)[which.max(best_r2)]
r2_bestoverall_tune <- best_r2[which.max(best_r2)]


all_prs_tune <- as.data.frame(all_prs_tune)
all_prs_valid <- as.data.frame(all_prs_valid)

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

tune_dat_sl_R2 <- data.frame(y = y_tune,x = prs_best_tune)
valid_dat_sl_R2 <- data.frame(y = y_valid,x = prs_best_vad)

r2_sl_tune <- summary(lm(y~x,data = tune_dat_sl_R2))$r.squared

if(r2_sl_tune < r2_bestoverall_tune){
  prs_best_tune_sl <- all_tune[,best_thresh]
  prs_best_vad_sl <- all_valid[,best_thresh] 
  
  tune_dat_sl_R2 <- data.frame(y = y_tune,x = prs_best_tune_sl)
  valid_dat_sl_R2 <- data.frame(y = y_valid,x = prs_best_vad_sl)
}





if(r2_sl_tune < r2_bestoverall_tune){
  beta_matrix$Beta_Final <- beta_matrix[,best_thresh]
  write.csv(beta_matrix[,c("Gene","Chr","Category","Number_SNV","Beta_Final")],file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"),row.names = FALSE)
}else{
  if(best_algorithm == "SL.glmnet_All"){
    beta_full_superlearner <- data.frame(Coef = row.names(coef(full_superlearner$fitLibrary$SL.glmnet_All$object)),Beta = as.numeric(coef(full_superlearner$fitLibrary$SL.glmnet_All$object)))
  }else if(best_algorithm == "SL.ridge_All"){
    beta_full_superlearner <- data.frame(Coef = row.names(full_superlearner$fitLibrary$SL.ridge_All$bestCoef),Beta = as.numeric(full_superlearner$fitLibrary$SL.ridge_All$bestCoef[,1]))
  }else if(best_algorithm == "SL.glm_All"){
    beta_full_superlearner <- data.frame(Coef = names(coef(full_superlearner$fitLibrary$SL.glm_All$object)),Beta = as.numeric(coef(full_superlearner$fitLibrary$SL.glm_All$object)))
  }else{
    beta_lasso <- data.frame(Coef = row.names(coef(full_superlearner$fitLibrary$SL.glmnet_All$object)),Beta = as.numeric(coef(full_superlearner$fitLibrary$SL.glmnet_All$object)))
    beta_ridge <- data.frame(Coef = row.names(full_superlearner$fitLibrary$SL.ridge_All$bestCoef),Beta = as.numeric(full_superlearner$fitLibrary$SL.ridge_All$bestCoef[,1]))
    beta_glm <- data.frame(Coef = names(coef(full_superlearner$fitLibrary$SL.glm_All$object)),Beta = as.numeric(coef(full_superlearner$fitLibrary$SL.glm_All$object)))
    beta_mean <- data.frame(Coef = names(coef(full_superlearner$fitLibrary$SL.glm_All$object)), Beta = 0)
    beta_mean$Beta[1] <- full_superlearner$fitLibrary$SL.mean_All$object
    
    beta_full_superlearner <- data.frame(Coef = names(coef(full_superlearner$fitLibrary$SL.glm_All$object)), Beta = 0)
    beta_full_superlearner$Beta <- as.numeric(full_superlearner$coef[names(full_superlearner$coef) == "SL.glmnet_All"]) * beta_lasso$Beta + 
      as.numeric(full_superlearner$coef[names(full_superlearner$coef) == "SL.ridge_All"]) * beta_ridge$Beta + 
      as.numeric(full_superlearner$coef[names(full_superlearner$coef) == "SL.glm_All"]) * beta_glm$Beta + 
      as.numeric(full_superlearner$coef[names(full_superlearner$coef) == "SL.mean_All"]) * beta_mean$Beta 
  }
  
  beta_full_superlearner <- beta_full_superlearner[-1,]
  
  beta_matrix$Beta_Final <- 0
  for(i in 1:nrow(beta_full_superlearner)){
    beta_matrix$Beta_Final <- beta_matrix[,beta_full_superlearner$Coef[i]]*beta_full_superlearner$Beta[i]
  }
  write.csv(beta_matrix[,c("Gene","Chr","Category","Number_SNV","Beta_Final")],file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_final_coef.csv"),row.names = FALSE)
}





load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

RV_PRS <- data.frame(IID = pheno_valid$IID,Y = valid_dat_sl_R2[,1],RV_PRS = valid_dat_sl_R2[,2],pheno_valid[,c("pc1","pc2","pc3","pc4","pc5")])

write.csv(RV_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"_BestPRS.csv"),row.names = FALSE)

RV_PRS_raw <- RV_PRS
RV_PRS_adjusted <- RV_PRS


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

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","SAS","AMR","AFR","EAS"), 
                         beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_SAS,best_beta_raw_RV_AMR,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                         se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_SAS,se_beta_raw_RV_AMR,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                         beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_AMR,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                         se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_AMR,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS))

write.csv(RV_PRS_Results,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_RareVariants_PRS/",trait,"Best_Betas.csv"),row.names = FALSE)
