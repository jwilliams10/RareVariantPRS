rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)

#Train
prs_train_CT <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/CT/prs_all_train.txt", sep="")
prs_train_LDPred2 <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_train_prs_all.txt")
prs_train_LASSOSum2 <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_train_prs_all.txt")

prs_train_all <- cbind(prs_train_LDPred2[,1:2],prs_train_CT[,-c(1,2)],prs_train_LDPred2[,-c(1,2,3,4)],prs_train_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_train_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_train_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_train_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_train_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_train_CT);rm(prs_train_LDPred2);rm(prs_train_LASSOSum2)

pheno_train <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt")
colnames(pheno_train) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
## Merge covariates and y for tuning with the prs_mat
pheno_train <- left_join(pheno_train,prs_train_all,by = "IID")

prs_train_all <- pheno_train[,-c(1:17)]

#Tune
prs_tune_CT <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/CT/prs_all_tune.txt", sep="")
prs_tune_LDPred2 <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_tune_prs_all.txt")
prs_tune_LASSOSum2 <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_tune_prs_all.txt")

prs_tune_all <- cbind(prs_tune_LDPred2[,1:2],prs_tune_CT[,-c(1,2)],prs_tune_LDPred2[,-c(1,2,3,4)],prs_tune_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_tune_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_tune_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_tune_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_tune_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_tune_CT);rm(prs_tune_LDPred2);rm(prs_tune_LASSOSum2)

pheno_tune <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Tune.txt")
colnames(pheno_tune) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
## Merge covariates and y for tuning with the prs_mat
pheno_tune <- left_join(pheno_tune,prs_tune_all,by = "IID")

prs_tune_all <- pheno_tune[,-c(1:17)]

#Validation
prs_validation_CT <- read.csv("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/CT/prs_all_validation.txt", sep="")
prs_validation_LDPred2 <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LDPred2/ldpred2_validation_prs_all.txt")
prs_validation_LASSOSum2 <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/LASSOSUM2/LASSOSUM2_validation_prs_all.txt")

prs_validation_all <- cbind(prs_validation_LDPred2[,1:2],prs_validation_CT[,-c(1,2)],prs_validation_LDPred2[,-c(1,2,3,4)],prs_validation_LASSOSum2[,-c(1,2,3,4)])
colnames(prs_validation_all) <- c("X.FID","IID",paste0("CT_",colnames(prs_validation_CT[,-c(1,2)])),paste0("LDPred2_",colnames(prs_validation_LDPred2[,-c(1,2,3,4)])),paste0("LASSOSum2_",colnames(prs_validation_LASSOSum2[,-c(1,2,3,4)])))
rm(prs_validation_CT);rm(prs_validation_LDPred2);rm(prs_validation_LASSOSum2)

pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Validation.txt")
colnames(pheno_validation) <- c("IID","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
## Merge covariates and y for tuning with the prs_mat
pheno_validation <- left_join(pheno_validation,prs_validation_all,by = "IID")

prs_validation_all <- pheno_validation[,-c(1:17)]

## Null Models

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_train)
y_train <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_tune)
y_tune <- model.null$residual

model.null <- lm(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pheno_validation)
y_validation <- model.null$residual

##################################

mtx <- cor(prs_tune_all)
drop <- findCorrelation(mtx,cutoff=0.98)
drop <- names(prs_tune_all)[drop]

prs_train_all = prs_train_all %>% 
  select(-all_of(drop))
prs_tune_all = prs_tune_all %>% 
  select(-all_of(drop))
prs_validation_all = prs_validation_all %>% 
  select(-all_of(drop))


## SL

SL.library <- c(
  "SL.glmnet",
  "SL.ridge",
  "SL.glm",
  "SL.mean"
)
sl <- SuperLearner(Y = y_tune, X = prs_tune_all, family = gaussian(),
                   # For a real analysis we would use V = 10.
                   # V = 3,
                   SL.library = SL.library)
cvsl <- CV.SuperLearner(Y = y_tune, X = prs_tune_all, family = gaussian(),
                        # For a real analysis we would use V = 10.
                        # V = 3,
                        SL.library = SL.library)

best_algorithm <- summary(cvsl)$Table$Algorithm[which.min(summary(cvsl)$Table$Ave)]

### Extract Coefficients
#algorithm weight
alg_weights <- sl$coef
#glmnet
glmnet_obj <- sl$fitLibrary$SL.glmnet$object
best_lambda <- glmnet_obj$lambda[which.min(glmnet_obj$cvm)]
glmnet_coefs <- coef(glmnet_obj, s = best_lambda)
#ridge
ridge_coefs <- sl$fitLibrary$SL.ridge_All$bestCoef
#glm 
glm_coefs <- sl$fitLibrary$SL.glm_All$object$coefficients
#mean 
mean_coefs <- sl$fitLibrary$SL.mean_All$object

if(best_algorithm == "SL.glmnet_All"){
  #final
  final_coefs <- glmnet_coefs
}else if(best_algorithm == "SL.ridge_All"){
  #final
  final_coefs <- ridge_coefs
}else if(best_algorithm == "SL.glm_All"){
  #final
  final_coefs <- glm_coefs
}else{
  #final
  final_coefs <- alg_weights[1] * glmnet_coefs + alg_weights[2] * ridge_coefs + alg_weights[3] * glm_coefs + alg_weights[4] * mean_coefs 
}

#remove the intercept
final_coefs = final_coefs[2:nrow(final_coefs),]
#remove weight 0 coefficients
final_coefs = final_coefs[final_coefs!=0]

save(final_coefs,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/final_coefs_All.RData")

a <- predict(sl, prs_validation_all, onlySL = FALSE)

prs_best_validation_sl <- a$pred
prs_best_validation_glmnet <- a$library.predict[,1]
prs_best_validation_ridge <- a$library.predict[,2]
prs_best_validation_glm <- a$library.predict[,3]

if(best_algorithm == "SL.glmnet_All"){
  #final
  prs_best_validation <- prs_best_validation_glmnet
}else if(best_algorithm == "SL.ridge_All"){
  #final
  prs_best_validation <- prs_best_validation_ridge
}else if(best_algorithm == "SL.glm_All"){
  #final
  prs_best_validation <- prs_best_validation_glm
}else{
  #final
  prs_best_validation <- prs_best_validation_sl
}

model <- lm(y_validation~prs_best_validation)
r2 <- summary(model)$r.square

prs_best_validation <- data.frame(IID = pheno_validation$IID,prs = prs_best_validation)


a <- predict(sl, prs_train_all, onlySL = FALSE)

prs_best_train_sl <- a$pred
prs_best_train_glmnet <- a$library.predict[,1]
prs_best_train_ridge <- a$library.predict[,2]
prs_best_train_glm <- a$library.predict[,3]

if(best_algorithm == "SL.glmnet_All"){
  #final
  prs_best_train <- prs_best_train_glmnet
}else if(best_algorithm == "SL.ridge_All"){
  #final
  prs_best_train <- prs_best_train_ridge
}else if(best_algorithm == "SL.glm_All"){
  #final
  prs_best_train <- prs_best_train_glm
}else{
  #final
  prs_best_train <- prs_best_train_sl
}
prs_best_train <- data.frame(IID = pheno_train$IID,prs = prs_best_train)


a <- predict(sl, prs_tune_all, onlySL = FALSE)

prs_best_tune_sl <- a$pred
prs_best_tune_glmnet <- a$library.predict[,1]
prs_best_tune_ridge <- a$library.predict[,2]
prs_best_tune_glm <- a$library.predict[,3]

if(best_algorithm == "SL.glmnet_All"){
  #final
  prs_best_tune <- prs_best_tune_glmnet
}else if(best_algorithm == "SL.ridge_All"){
  #final
  prs_best_tune <- prs_best_tune_ridge
}else if(best_algorithm == "SL.glm_All"){
  #final
  prs_best_tune <- prs_best_tune_glm
}else{
  #final
  prs_best_tune <- prs_best_tune_sl
}
prs_best_tune <- data.frame(IID = pheno_tune$IID,prs = prs_best_tune)

write.table(prs_best_train,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Train_All.txt",sep = "\t",row.names = FALSE)
write.table(prs_best_tune,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Tune_All.txt",sep = "\t",row.names = FALSE)
write.table(prs_best_validation,file="/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/Best_Validation_All.txt",sep = "\t",row.names = FALSE)

## bootstrap the R2 to provide an approximate distribution 
data <- data.frame(y = y_validation, x = prs_best_validation$prs)
R2Boot <- function(data,indices){
  boot_data <- data[indices, ]
  model <- lm(y ~ x, data = boot_data)
  result <- summary(model)$r.square
  return(c(result))
}
boot_r2 <- boot(data = data, statistic = R2Boot, R = 1000)

ci_result <- boot.ci(boot_r2, type = "perc")
SL.result <- data.frame(method = "SL_Combined",
                        r2 = r2,
                        r2_low = ci_result$percent[4],
                        r2_high = ci_result$percent[5]
)

save(SL.result,file = "/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/Combined_Common_PRS/sl_result_All.RData")
