

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

system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --pheno /data/williamsjacr/UKB_WES_Phenotypes/All_Train_NewPCs.txt --pheno-name ",trait," --linear --covar /data/williamsjacr/UKB_WES_Phenotypes/All_Train_NewPCs.txt --covar-name age, age2, sex, ",paste(paste0("pc",1:10),collapse = ",")," --vif 999 --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats_newpcs"))

