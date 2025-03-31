rm(list = ls())

library(bigsnpr)

system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/reference.txt --maf 0.01 --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference"))

if(file.exists("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.bk")){
  file.remove("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.bk")
}
if(file.exists("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.rds")){
  file.remove("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.rds")
}
snp_readBed("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.bed",backingfile = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference")