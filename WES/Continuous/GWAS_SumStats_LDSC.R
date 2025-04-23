

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

## plink 10 PCs
system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --pheno /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --pheno-name ",trait," --linear --covar /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --covar-name age, age2, sex, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10 --vif 999 --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats"))

## plink 20 PCs
system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --pheno /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --pheno-name ",trait," --linear --covar /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --covar-name age, age2, sex, ",paste(paste0("pc",1:20),collapse = ",")," --vif 999 --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats_20PCs"))

## plink 40 PCs
system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --pheno /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --pheno-name ",trait," --linear --covar /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --covar-name age, age2, sex, ",paste(paste0("pc",1:40),collapse = ",")," --vif 999 --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats_40PCs"))

## plink ranknormal
system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --pheno /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --pheno-name ",trait,"adj_norm --linear --covar /data/williamsjacr/UKB_WES_Phenotypes/All_Train_LDSC.txt --covar-name age, age2, sex, ",paste(paste0("pc",1:10),collapse = ",")," --vif 999 --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats_ranknormal"))
