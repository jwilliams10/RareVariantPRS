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

Gene_Centric_Coding_G_Star <- function(chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous"),
                                       genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                       QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                       Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){
  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  genes <- genes_info[genes_info[,2]==chr,]  
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  ## get SNV id, position, REF, ALT (whole genome)
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant"){
    SNVlist <- filter == "PASS"
  }
  
  if(variant_type=="SNV"){
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  
  if(variant_type=="Indel"){
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  
  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  
  rm(filter)
  gc()
  
  ### Gene
  kk <- which(genes[,1]==gene_name)
  
  sub_start_loc <- genes[kk,3]
  sub_end_loc <- genes[kk,4]
  
  is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
  variant.id.gene <- variant.id[is.in]
  
  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
  
  ## plof
  ## Gencode_Exonic
  GENCODE.EXONIC.Category  <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
  
  variant.id.gene <- seqGetData(genofile, "variant.id")  
  
  if(category == "plof"){
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
    variant.id.gene <- variant.id.gene[lof.in.plof]
  }else if(category == "plof_ds"){
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene <- variant.id.gene[lof.in.plof]
  }else if(category == "missense"){
    lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.missense]
  }else if(category == "disruptive_missense"){
    lof.in.ds <- ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene <- variant.id.gene[lof.in.ds]
  }else{
    lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.synonymous]
  }
  
  seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ## Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno)))
  {
    if(dim(Geno)[2]>0)
    {
      if(geno_missing_imputation=="mean")
      {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor")
      {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  genotype <- Geno
  
  if(dim(genotype)[2] == 1){
    return(matrix(0,nrow = dim(genotype)[1],ncol = 1))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]
  G <- Geno_rare
  rm(Geno_rare)
  gc()
  
  if(is.null(dim(G))){
    G <- matrix(G,ncol = 1)
  }
  
  C <- G%*%matrix(1,nrow=ncol(G),ncol = 1)
  
  seqResetFilter(genofile)
  
  return(C)
}



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

Train_PVals_All <- read.csv(paste0(trait,"_coding_sig.csv"))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= 1e-03,]

## agds dir
agds_dir <- paste0("ukb.200k.wgs.chr",1:22,".pass.annotated.gds")

## Null Model
obj_nullmodel_train <- get(load(paste0(trait,"_Train_Null_Model.RData")))
obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)

## Parameter
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("Annotation_name_catalog.Rdata"))

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
pheno_train <- read.delim("All_Train.txt")
pheno_train <- inner_join(pheno_train,X_train)
X_train <- as.matrix(pheno_train[,27:ncol(pheno_train),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_train)
y_train <- model.null$residual

X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_tune$id_include],G_star_gene_centric_coding_tune)
pheno_tune <- read.delim("All_Tune.txt")
pheno_tune <- inner_join(pheno_tune,X_tune)
X_tune <- as.matrix(pheno_tune[,27:ncol(pheno_tune),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
y_tune <- model.null$residual

X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_coding_vad)
pheno_valid <- read.delim("All_Validation.txt")
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
  write.csv(beta_matrix[,c("Gene","Chr","Category","Number_SNV","Beta_Final")],file = paste0(trait,"_Final_Coefficients_Coding.csv"),row.names = FALSE)
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
  write.csv(beta_matrix[,c("Gene","Chr","Category","Number_SNV","Beta_Final")],file = paste0(trait,"_Final_Coefficients_Coding.csv"),row.names = FALSE)
}





load("all_phenotypes.RData")

RV_PRS_Coding <- data.frame(IID = pheno_valid$IID,Y = valid_dat_sl_R2[,1],RV_PRS = valid_dat_sl_R2[,2],pheno_valid[,c("pc1","pc2","pc3","pc4","pc5")])

write.csv(RV_PRS,file = paste0(trait,"_BestPRS_Coding.csv"),row.names = FALSE)

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
RV_PRS_raw_NonEUR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
RV_PRS_raw_UNK <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
RV_PRS_raw_SAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
RV_PRS_raw_MIX <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
RV_PRS_raw_AFR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
RV_PRS_raw_EAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

RV_PRS_adjusted_EUR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
RV_PRS_adjusted_NonEUR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
RV_PRS_adjusted_UNK <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
RV_PRS_adjusted_SAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
RV_PRS_adjusted_MIX <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
RV_PRS_adjusted_AFR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
RV_PRS_adjusted_EAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

RV_PRS_raw_EUR$Y <- scale(RV_PRS_raw_EUR$Y)
RV_PRS_raw_NonEUR$Y <- scale(RV_PRS_raw_NonEUR$Y)
RV_PRS_raw_UNK$Y <- scale(RV_PRS_raw_UNK$Y)
RV_PRS_raw_SAS$Y <- scale(RV_PRS_raw_SAS$Y)
RV_PRS_raw_MIX$Y <- scale(RV_PRS_raw_MIX$Y)
RV_PRS_raw_AFR$Y <- scale(RV_PRS_raw_AFR$Y)
RV_PRS_raw_EAS$Y <- scale(RV_PRS_raw_EAS$Y)

RV_PRS_adjusted_EUR$Y <- scale(RV_PRS_adjusted_EUR$Y)
RV_PRS_adjusted_NonEUR$Y <- scale(RV_PRS_adjusted_NonEUR$Y)
RV_PRS_adjusted_UNK$Y <- scale(RV_PRS_adjusted_UNK$Y)
RV_PRS_adjusted_SAS$Y <- scale(RV_PRS_adjusted_SAS$Y)
RV_PRS_adjusted_MIX$Y <- scale(RV_PRS_adjusted_MIX$Y)
RV_PRS_adjusted_AFR$Y <- scale(RV_PRS_adjusted_AFR$Y)
RV_PRS_adjusted_EAS$Y <- scale(RV_PRS_adjusted_EAS$Y)

RV_PRS_raw_EUR$RV_PRS <- scale(RV_PRS_raw_EUR$RV_PRS)
RV_PRS_raw_NonEUR$RV_PRS <- scale(RV_PRS_raw_NonEUR$RV_PRS)
RV_PRS_raw_UNK$RV_PRS <- scale(RV_PRS_raw_UNK$RV_PRS)
RV_PRS_raw_SAS$RV_PRS <- scale(RV_PRS_raw_SAS$RV_PRS)
RV_PRS_raw_MIX$RV_PRS <- scale(RV_PRS_raw_MIX$RV_PRS)
RV_PRS_raw_AFR$RV_PRS <- scale(RV_PRS_raw_AFR$RV_PRS)
RV_PRS_raw_EAS$RV_PRS <- scale(RV_PRS_raw_EAS$RV_PRS)


best_beta_raw_RV_EUR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_EUR))[2]
se_beta_raw_RV_EUR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_EUR))$coefficients[2,2]
best_beta_raw_RV_NonEUR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_NonEUR))[2]
se_beta_raw_RV_NonEUR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_NonEUR))$coefficients[2,2]
best_beta_raw_RV_UNK <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_UNK))[2]
se_beta_raw_RV_UNK <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_UNK))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_SAS))[2]
se_beta_raw_RV_SAS <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_SAS))$coefficients[2,2]
best_beta_raw_RV_MIX <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_MIX))[2]
se_beta_raw_RV_MIX <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_MIX))$coefficients[2,2]
best_beta_raw_RV_AFR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_AFR))[2]
se_beta_raw_RV_AFR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_AFR))$coefficients[2,2]
best_beta_raw_RV_EAS <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_EAS))[2]
se_beta_raw_RV_EAS <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_EAS))$coefficients[2,2]

best_beta_adjusted_RV_EUR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_EUR))[2]
se_beta_adjusted_RV_EUR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_EUR))$coefficients[2,2]
best_beta_adjusted_RV_NonEUR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_NonEUR))[2]
se_beta_adjusted_RV_NonEUR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_NonEUR))$coefficients[2,2]
best_beta_adjusted_RV_UNK <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_UNK))[2]
se_beta_adjusted_RV_UNK <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_UNK))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_SAS))[2]
se_beta_adjusted_RV_SAS <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_SAS))$coefficients[2,2]
best_beta_adjusted_RV_MIX <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_MIX))[2]
se_beta_adjusted_RV_MIX <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_MIX))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_AFR))[2]
se_beta_adjusted_RV_AFR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_AFR))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_EAS))[2]
se_beta_adjusted_RV_EAS <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_EAS))$coefficients[2,2]

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_NonEUR,best_beta_raw_RV_UNK,best_beta_raw_RV_SAS,best_beta_raw_RV_MIX,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_NonEUR,se_beta_raw_RV_UNK,se_beta_raw_RV_SAS,se_beta_raw_RV_MIX,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_NonEUR,best_beta_adjusted_RV_UNK,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_MIX,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                             se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_NonEUR,se_beta_adjusted_RV_UNK,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_MIX,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS))

write.csv(RV_PRS_Results,file = paste0(trait,"Best_Betas_Coding.csv"),row.names = FALSE)

rm(list=setdiff(ls(), "trait"))































Gene_Centric_Noncoding_G_Star <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS","ncRNA"),
                                          genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                          QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                          Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,silent=FALSE){
  
  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  
  genes <- genes_info[genes_info[,2]==chr,]
  
  phenotype.id <- as.character(obj_nullmodel$id_include)
  
  if(category == "downstream"){
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    
    rm(filter)
    
    ## downstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="downstream")&(SNVlist)
    variant.id.downstream <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)
    rm(variant.id.downstream)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
    
    rm(GENCODE.Info)
    rm(variant_gene_num)
    
    Gene <- as.character(unlist(GENCODE.Info.split))
    
    rm(GENCODE.Info.split)
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category == "upstream"){
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    
    rm(filter)
    
    ## upstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- (GENCODE.Category=="upstream")&(SNVlist)
    variant.id.upstream <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)
    rm(variant.id.upstream)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
    variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
    
    rm(GENCODE.Info)
    rm(variant_gene_num)
    
    Gene <- as.character(unlist(GENCODE.Info.split))
    rm(GENCODE.Info.split)
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category == "UTR"){
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    rm(filter)
    
    ## downstream SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
    variant.id.UTR <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)
    rm(variant.id.UTR)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
    
    rm(GENCODE.Info)
    
    # Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[seq(1,length(z),2)]))
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))
    rm(GENCODE.Info.split)
    
    variant.id.SNV <- seqGetData(genofile, "variant.id")
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- which(Gene==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category == "promoter_CAGE"){
    ## Promoter
    varid <- seqGetData(genofile, "variant.id")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)
    
    #Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGEBvt <- CAGEAnno!=""
    CAGEidx <- which(CAGEBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEidx])
    seqSetFilter(genofile,promGobj,intersect=TRUE)
    CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
    ##obtain variants info
    CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
    CAGEvref <- as.character(seqGetData(genofile,"$ref"))
    CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)
    
    rm(varid)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
    dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
    dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
    dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)
    
    seqResetFilter(genofile)
    rm(dfPromCAGEVarGene)
    
    ### Gene
    is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category=="promoter_DHS"){
    ## Promoter
    varid <- seqGetData(genofile, "variant.id")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)
    
    # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRsBvt <- rOCRsAnno!=""
    rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsidx])
    
    seqSetFilter(genofile,promGobj,intersect=TRUE)
    rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
    ## obtain variants info
    rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
    rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
    rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
    rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
    dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)
    
    rm(varid)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
    dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
    dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
    dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)
    
    seqResetFilter(genofile)
    rm(dfPromrOCRsVarGene)
    
    ### Gene
    is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category=="enhancer_CAGE"){
    ## Enhancer
    varid <- seqGetData(genofile, "variant.id")
    
    #Now extract the GeneHancer with CAGE Signal Overlay
    genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    genehancer <- genehancerAnno!=""
    
    CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
    CAGE <- CAGEAnno!=""
    CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
    CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])
    
    # variants that covered by whole GeneHancer without CAGE overlap.
    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerCAGEVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)
    
    rm(varid)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist,]
    dfHancerCAGEVarGene.SNV$enhancervpos <- as.character(dfHancerCAGEVarGene.SNV$enhancervpos)
    dfHancerCAGEVarGene.SNV$enhancervref <- as.character(dfHancerCAGEVarGene.SNV$enhancervref)
    dfHancerCAGEVarGene.SNV$enhancervalt <- as.character(dfHancerCAGEVarGene.SNV$enhancervalt)
    
    seqResetFilter(genofile)
    
    rm(dfHancerCAGEVarGene)
    
    ### Gene
    is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else if(category=="enhancer_DHS"){
    ## Enhancer
    varid <- seqGetData(genofile, "variant.id")
    
    #Now extract the GeneHancer with rOCRs Signal Overlay
    genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    genehancer <- genehancerAnno!=""
    
    rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
    rOCRs <- rOCRsAnno!=""
    rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
    rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
    seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
    # variants that covered by whole GeneHancer without rOCRs overlap.
    
    genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
    enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
    enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
    enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
    enhancervpos <- as.numeric(seqGetData(genofile,"position"))
    enhancervref <- as.character(seqGetData(genofile,"$ref"))
    enhancervalt <- as.character(seqGetData(genofile,"$alt"))
    dfHancerrOCRsVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)
    
    rm(varid)
    
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    variant.id.SNV <- variant.id[SNVlist]
    
    dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist,]
    dfHancerrOCRsVarGene.SNV$enhancervpos <- as.character(dfHancerrOCRsVarGene.SNV$enhancervpos)
    dfHancerrOCRsVarGene.SNV$enhancervref <- as.character(dfHancerrOCRsVarGene.SNV$enhancervref)
    dfHancerrOCRsVarGene.SNV$enhancervalt <- as.character(dfHancerrOCRsVarGene.SNV$enhancervalt)
    
    seqResetFilter(genofile)
    rm(dfHancerrOCRsVarGene)
    
    ### Gene
    is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
    variant.is.in <- variant.id.SNV[is.in]
    
  }else{
    ## get SNV id
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant"){
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV"){
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel"){
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    variant.id <- seqGetData(genofile, "variant.id")
    rm(filter)
    
    ## ncRNA SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)
    
    variant.id.ncRNA <- variant.id[is.in]
    rm(GENCODE.Category)
    
    seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)
    rm(variant.id.ncRNA)
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))
    
    Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
    Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
    Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))
    
    rm(GENCODE.Info)
    rm(GENCODE.Info.split)
    
    variant.id.ncRNA <- seqGetData(genofile, "variant.id")
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
    is.in <- union(is.in,which(Gene_list_3==gene_name))
    
    variant.is.in <- variant.id.ncRNA[is.in]
    
  }
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))
  
  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  ##Genotype
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match,,drop=FALSE]
  
  ## impute missing
  if(!is.null(dim(Geno))){
    if(dim(Geno)[2]>0){
      if(geno_missing_imputation=="mean"){
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if(geno_missing_imputation=="minor"){
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  
  genotype <- Geno
  
  if(dim(genotype)[2] == 1){
    return(matrix(0,nrow = dim(genotype)[1],ncol = 1))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]
  G <- Geno_rare
  rm(Geno_rare)
  gc()
  
  if(is.null(dim(G))){
    G <- matrix(G,ncol = 1)
  }
  
  C <- G%*%matrix(1,nrow=ncol(G),ncol = 1)
  
  seqResetFilter(genofile)
  
  return(C)
}


Train_PVals_All <- read.csv(paste0(trait,"_noncoding_sig.csv"))
Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= 1e-03,]

## agds dir
agds_dir <- paste0("ukb.200k.wgs.chr",1:22,".pass.annotated.gds")

## Null Model
obj_nullmodel_train <- get(load(paste0(trait,"_Train_Null_Model.RData")))
obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)

## Parameter
QC_label <- "annotation/info/QC_label2"
geno_missing_imputation <- "mean"
variant_type <- "SNV"	

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")

G_star_gene_centric_noncoding <- list()

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
  
  G_star_gene_centric_noncoding[[i]] <- Gene_Centric_Noncoding_G_Star(chr=chr,gene_name=gene_name,category=category ,
                                                                genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent=FALSE) 
  seqClose(genofile) 
} 

G_star_gene_centric_noncoding <- do.call(cbind,G_star_gene_centric_noncoding)

col_remove <- apply(G_star_gene_centric_noncoding,2,function(x){sum(x != 0)}) > 10 & colSums(G_star_gene_centric_noncoding) > 10 
G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[,col_remove,drop = FALSE]

Train_PVals_All <- Train_PVals_All[col_remove,]

ids_gstar <- obj_nullmodel$id_include

G_star_gene_centric_noncoding_train <- G_star_gene_centric_noncoding[ids_gstar %in% obj_nullmodel_train$id_include,]
G_star_gene_centric_noncoding_tune <- G_star_gene_centric_noncoding[ids_gstar %in% obj_nullmodel_tune$id_include,]
G_star_gene_centric_noncoding_vad <- G_star_gene_centric_noncoding[ids_gstar %in% obj_nullmodel_validation$id_include,]

rm(G_star_gene_centric_noncoding)

X_train <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_train$id_include],G_star_gene_centric_noncoding_train)
pheno_train <- read.delim("All_Train.txt")
pheno_train <- inner_join(pheno_train,X_train)
X_train <- as.matrix(pheno_train[,27:ncol(pheno_train),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_train)
y_train <- model.null$residual

X_tune <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_tune$id_include],G_star_gene_centric_noncoding_tune)
pheno_tune <- read.delim("All_Tune.txt")
pheno_tune <- inner_join(pheno_tune,X_tune)
X_tune <- as.matrix(pheno_tune[,27:ncol(pheno_tune),drop = FALSE])
model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_tune)
y_tune <- model.null$residual

X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_noncoding_vad)
pheno_valid <- read.delim("All_Validation.txt")
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
  write.csv(beta_matrix[,c("Gene","Chr","Category","Number_SNV","Beta_Final")],file = paste0(trait,"_Final_Coefficients_Coding.csv"),row.names = FALSE)
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
  write.csv(beta_matrix[,c("Gene","Chr","Category","Number_SNV","Beta_Final")],file = paste0(trait,"_Final_Coefficients_Noncoding.csv"),row.names = FALSE)
}





load("all_phenotypes.RData")

RV_PRS_Coding <- data.frame(IID = pheno_valid$IID,Y = valid_dat_sl_R2[,1],RV_PRS = valid_dat_sl_R2[,2],pheno_valid[,c("pc1","pc2","pc3","pc4","pc5")])

write.csv(RV_PRS,file = paste0(trait,"_BestPRS_Noncoding.csv"),row.names = FALSE)

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
RV_PRS_raw_NonEUR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
RV_PRS_raw_UNK <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
RV_PRS_raw_SAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
RV_PRS_raw_MIX <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
RV_PRS_raw_AFR <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
RV_PRS_raw_EAS <- RV_PRS_raw[RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

RV_PRS_adjusted_EUR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
RV_PRS_adjusted_NonEUR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
RV_PRS_adjusted_UNK <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
RV_PRS_adjusted_SAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
RV_PRS_adjusted_MIX <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
RV_PRS_adjusted_AFR <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
RV_PRS_adjusted_EAS <- RV_PRS_adjusted[RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

RV_PRS_raw_EUR$Y <- scale(RV_PRS_raw_EUR$Y)
RV_PRS_raw_NonEUR$Y <- scale(RV_PRS_raw_NonEUR$Y)
RV_PRS_raw_UNK$Y <- scale(RV_PRS_raw_UNK$Y)
RV_PRS_raw_SAS$Y <- scale(RV_PRS_raw_SAS$Y)
RV_PRS_raw_MIX$Y <- scale(RV_PRS_raw_MIX$Y)
RV_PRS_raw_AFR$Y <- scale(RV_PRS_raw_AFR$Y)
RV_PRS_raw_EAS$Y <- scale(RV_PRS_raw_EAS$Y)

RV_PRS_adjusted_EUR$Y <- scale(RV_PRS_adjusted_EUR$Y)
RV_PRS_adjusted_NonEUR$Y <- scale(RV_PRS_adjusted_NonEUR$Y)
RV_PRS_adjusted_UNK$Y <- scale(RV_PRS_adjusted_UNK$Y)
RV_PRS_adjusted_SAS$Y <- scale(RV_PRS_adjusted_SAS$Y)
RV_PRS_adjusted_MIX$Y <- scale(RV_PRS_adjusted_MIX$Y)
RV_PRS_adjusted_AFR$Y <- scale(RV_PRS_adjusted_AFR$Y)
RV_PRS_adjusted_EAS$Y <- scale(RV_PRS_adjusted_EAS$Y)

RV_PRS_raw_EUR$RV_PRS <- scale(RV_PRS_raw_EUR$RV_PRS)
RV_PRS_raw_NonEUR$RV_PRS <- scale(RV_PRS_raw_NonEUR$RV_PRS)
RV_PRS_raw_UNK$RV_PRS <- scale(RV_PRS_raw_UNK$RV_PRS)
RV_PRS_raw_SAS$RV_PRS <- scale(RV_PRS_raw_SAS$RV_PRS)
RV_PRS_raw_MIX$RV_PRS <- scale(RV_PRS_raw_MIX$RV_PRS)
RV_PRS_raw_AFR$RV_PRS <- scale(RV_PRS_raw_AFR$RV_PRS)
RV_PRS_raw_EAS$RV_PRS <- scale(RV_PRS_raw_EAS$RV_PRS)


best_beta_raw_RV_EUR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_EUR))[2]
se_beta_raw_RV_EUR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_EUR))$coefficients[2,2]
best_beta_raw_RV_NonEUR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_NonEUR))[2]
se_beta_raw_RV_NonEUR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_NonEUR))$coefficients[2,2]
best_beta_raw_RV_UNK <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_UNK))[2]
se_beta_raw_RV_UNK <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_UNK))$coefficients[2,2]
best_beta_raw_RV_SAS <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_SAS))[2]
se_beta_raw_RV_SAS <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_SAS))$coefficients[2,2]
best_beta_raw_RV_MIX <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_MIX))[2]
se_beta_raw_RV_MIX <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_MIX))$coefficients[2,2]
best_beta_raw_RV_AFR <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_AFR))[2]
se_beta_raw_RV_AFR <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_AFR))$coefficients[2,2]
best_beta_raw_RV_EAS <- coef(lm(Y~RV_PRS,data = RV_PRS_raw_EAS))[2]
se_beta_raw_RV_EAS <- summary(lm(Y~RV_PRS,data = RV_PRS_raw_EAS))$coefficients[2,2]

best_beta_adjusted_RV_EUR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_EUR))[2]
se_beta_adjusted_RV_EUR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_EUR))$coefficients[2,2]
best_beta_adjusted_RV_NonEUR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_NonEUR))[2]
se_beta_adjusted_RV_NonEUR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_NonEUR))$coefficients[2,2]
best_beta_adjusted_RV_UNK <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_UNK))[2]
se_beta_adjusted_RV_UNK <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_UNK))$coefficients[2,2]
best_beta_adjusted_RV_SAS <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_SAS))[2]
se_beta_adjusted_RV_SAS <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_SAS))$coefficients[2,2]
best_beta_adjusted_RV_MIX <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_MIX))[2]
se_beta_adjusted_RV_MIX <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_MIX))$coefficients[2,2]
best_beta_adjusted_RV_AFR <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_AFR))[2]
se_beta_adjusted_RV_AFR <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_AFR))$coefficients[2,2]
best_beta_adjusted_RV_EAS <- coef(lm(Y~RV_PRS,data = RV_PRS_adjusted_EAS))[2]
se_beta_adjusted_RV_EAS <- summary(lm(Y~RV_PRS,data = RV_PRS_adjusted_EAS))$coefficients[2,2]

RV_PRS_Results <- data.frame(trait = trait,ancestry = c("EUR","NonEUR","UNK","SAS","MIX","AFR","EAS"), 
                             beta_raw = c(best_beta_raw_RV_EUR,best_beta_raw_RV_NonEUR,best_beta_raw_RV_UNK,best_beta_raw_RV_SAS,best_beta_raw_RV_MIX,best_beta_raw_RV_AFR,best_beta_raw_RV_EAS), 
                             se_raw = c(se_beta_raw_RV_EUR,se_beta_raw_RV_NonEUR,se_beta_raw_RV_UNK,se_beta_raw_RV_SAS,se_beta_raw_RV_MIX,se_beta_raw_RV_AFR,se_beta_raw_RV_EAS), 
                             beta_adjusted = c(best_beta_adjusted_RV_EUR,best_beta_adjusted_RV_NonEUR,best_beta_adjusted_RV_UNK,best_beta_adjusted_RV_SAS,best_beta_adjusted_RV_MIX,best_beta_adjusted_RV_AFR,best_beta_adjusted_RV_EAS), 
                             se_adjusted = c(se_beta_adjusted_RV_EUR,se_beta_adjusted_RV_NonEUR,se_beta_adjusted_RV_UNK,se_beta_adjusted_RV_SAS,se_beta_adjusted_RV_MIX,se_beta_adjusted_RV_AFR,se_beta_adjusted_RV_EAS))

write.csv(RV_PRS_Results,file = paste0(trait,"Best_Betas_Noncoding.csv"),row.names = FALSE)
