rm(list = ls())
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

# for array in {1..5};
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/RareVariantCoefficients_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/RareVariantCoefficients_Binary.sh -icmd="bash RareVariantCoefficients_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/ --priority low --instance-type mem3_ssd1_v2_x4
# done

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

Coding_Train_PVals_All <- read.csv("coding_sig.csv")
system("rm coding_sig.csv")
Coding_Train_PVals_All <- Coding_Train_PVals_All[Coding_Train_PVals_All$Trait == trait,]
Coding_Train_PVals_All <- subset(Coding_Train_PVals_All,select = -Trait)
Coding_Train_PVals_All <- Coding_Train_PVals_All[Coding_Train_PVals_All$STAARB <= 1e-03,]

## agds dir

## Null Model
obj_nullmodel_train <- get(load(paste0(trait,"_Train_Null_Model.RData")))
obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)

chrs <- unique(Coding_Train_PVals_All$Chr)
G_star_gene_centric_coding <- read.csv(paste0(trait,"_G_Star_Coding_Chr",chrs[1],".csv"))
system(paste0("rm ",paste0(trait,"_G_Star_Coding_Chr",chrs[1],".csv")))
for(i in 2:length(chrs)){
  G_star_gene_centric_coding <- cbind(G_star_gene_centric_coding,read.csv(paste0(trait,"_G_Star_Coding_Chr",chrs[i],".csv")))
  system(paste0("rm ",paste0(trait,"_G_Star_Coding_Chr",chrs[i],".csv")))
}

col_remove <- apply(G_star_gene_centric_coding,2,function(x){sum(x != 0)}) > 10 & colSums(G_star_gene_centric_coding) > 10 
G_star_gene_centric_coding <- G_star_gene_centric_coding[,col_remove,drop = FALSE]

Coding_Train_PVals_All <- Coding_Train_PVals_All[col_remove,]

ids_gstar <- obj_nullmodel$id_include

G_star_gene_centric_coding_train <- G_star_gene_centric_coding[ids_gstar %in% obj_nullmodel_train$id_include,]
G_star_gene_centric_coding_tune <- G_star_gene_centric_coding[ids_gstar %in% obj_nullmodel_tune$id_include,]
G_star_gene_centric_coding_vad <- G_star_gene_centric_coding[ids_gstar %in% obj_nullmodel_validation$id_include,]

rm(G_star_gene_centric_coding)





Noncoding_Train_PVals_All <- read.csv("noncoding_sig.csv")
system("rm noncoding_sig.csv")
Noncoding_Train_PVals_All <- Noncoding_Train_PVals_All[Noncoding_Train_PVals_All$Trait == trait,]
Noncoding_Train_PVals_All <- subset(Noncoding_Train_PVals_All,select = -Trait)
Noncoding_Train_PVals_All <- Noncoding_Train_PVals_All[Noncoding_Train_PVals_All$STAARB <= 1e-03,]

## agds dir

## Null Model
obj_nullmodel_train <- get(load(paste0(trait,"_Train_Null_Model.RData")))
system(paste0("rm ",paste0(trait,"_Train_Null_Model.RData")))
obj_nullmodel_tune <- get(load(paste0(trait,"_Tune_Null_Model.RData")))
system(paste0("rm ",paste0(trait,"_Tune_Null_Model.RData")))
obj_nullmodel_validation <- get(load(paste0(trait,"_Validation_Null_Model.RData")))
system(paste0("rm ",paste0(trait,"_Validation_Null_Model.RData")))

obj_nullmodel <- obj_nullmodel_train
obj_nullmodel$id_include <- c(obj_nullmodel_train$id_include,obj_nullmodel_tune$id_include,obj_nullmodel_validation$id_include)


chrs <- unique(Noncoding_Train_PVals_All$Chr)
G_star_gene_centric_noncoding <- read.csv(paste0(trait,"_G_Star_Noncoding_Chr",chrs[1],".csv"))
system(paste0("rm ",paste0(trait,"_G_Star_Noncoding_Chr",chrs[1],".csv")))
for(i in 2:length(chrs)){
  G_star_gene_centric_noncoding <- cbind(G_star_gene_centric_noncoding,read.csv(paste0(trait,"_G_Star_Noncoding_Chr",chrs[i],".csv")))
  system(paste0("rm ",paste0(trait,"_G_Star_Noncoding_Chr",chrs[i],".csv")))
}

col_remove <- apply(G_star_gene_centric_noncoding,2,function(x){sum(x != 0)}) > 10 & colSums(G_star_gene_centric_noncoding) > 10 
G_star_gene_centric_noncoding <- G_star_gene_centric_noncoding[,col_remove,drop = FALSE]

Noncoding_Train_PVals_All <- Noncoding_Train_PVals_All[col_remove,]

ids_gstar <- obj_nullmodel$id_include

G_star_gene_centric_noncoding_train <- G_star_gene_centric_noncoding[ids_gstar %in% obj_nullmodel_train$id_include,]
G_star_gene_centric_noncoding_tune <- G_star_gene_centric_noncoding[ids_gstar %in% obj_nullmodel_tune$id_include,]
G_star_gene_centric_noncoding_vad <- G_star_gene_centric_noncoding[ids_gstar %in% obj_nullmodel_validation$id_include,]


Train_PVals_All <- rbind(Coding_Train_PVals_All,Noncoding_Train_PVals_All)

X_valid <- data.frame(IID = ids_gstar[ids_gstar %in% obj_nullmodel_validation$id_include],G_star_gene_centric_coding_vad,G_star_gene_centric_noncoding_vad)
RV_PRS <- read.csv(paste0(trait,"_BestPRS.csv"))
tmp <- inner_join(RV_PRS[,c("IID","RV_PRS")],X_valid)
tmp <- subset(tmp,select = -c(IID))
print(summary(lm(RV_PRS~.,tmp))$r.squared)
Train_PVals_All$Beta <- coef(lm(RV_PRS~.,tmp))[-1]
Train_PVals_All$Beta[is.na(Train_PVals_All$Beta)] <- 0
Train_PVals_All$Beta[abs(Train_PVals_All$Beta) < 1e-10] <- 0

write.csv(Train_PVals_All,file = paste0(trait,"_final_coef.csv"),row.names = FALSE)

system("rm All_Train.txt")
system("rm All_Tune.txt")
system("rm All_Validation.txt")
system("rm all_phenotypes.RData")
system(paste0("rm ",paste0(trait,"_BestPRS.csv")))
