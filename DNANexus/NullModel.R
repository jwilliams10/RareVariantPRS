rm(list = ls())

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)
library(SCANG)
library(STAARpipeline)


# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/NullModel.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/NullModel.sh -icmd="bash NullModel.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/ --priority low --instance-type mem1_ssd1_v2_x4

trait <- "BMI"

for(trait in c("BMI","TC","LDL","HDL","logTG","Height")){
  pheno_train <- read.delim("All_Train.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Train_All.txt"))
  
  pheno_train <- inner_join(pheno_train,common_prs,by = "IID")
  
  obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + prs")), data = pheno_train,id = "IID",kins = NULL,family = gaussian(link = "identity"))
  
  save(obj.STAAR.UKB,file = paste0(trait,"_Train_Null_Model.RData"))
  
  
  
  
  pheno_tune <- read.delim("All_Tune.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Tune_All.txt"))
  
  pheno_tune <- inner_join(pheno_tune,common_prs,by = "IID")
  
  obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + prs")), data = pheno_tune,id = "IID",kins = NULL,family = gaussian(link = "identity"))
  
  save(obj.STAAR.UKB,file = paste0(trait,"_Tune_Null_Model.RData"))
  
  
  
  
  pheno_validation <- read.delim("All_Validation.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Validation_All.txt"))
  
  pheno_validation <- inner_join(pheno_validation,common_prs,by = "IID")
  
  obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + prs")), data = pheno_validation,id = "IID",kins = NULL,family = gaussian(link = "identity"))
  
  save(obj.STAAR.UKB,file = paste0(trait,"_Validation_Null_Model.RData"))
  
  
  file.remove(paste0(trait,"_Best_Train_All.txt")) 
  file.remove(paste0(trait,"_Best_Tune_All.txt")) 
  file.remove(paste0(trait,"_Best_Validation_All.txt")) 
  
}

file.remove("All_Train.txt")
file.remove("All_Tune.txt")
file.remove("All_Validation.txt")