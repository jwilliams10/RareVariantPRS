

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --pheno /data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt --pheno-name ",trait," --linear --covar /data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt --covar-name age, age2, sex, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10 --vif 999 --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats"))
}
