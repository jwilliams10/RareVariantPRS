rm(list = ls())

library(data.table)
library(dplyr)

for(trait in c("Asthma","CAD","T2D","Breast","Prostate","BMI","LDL","HDL","logTG","TC","Height")){
  if(trait %in% c("BMI","LDL","HDL","logTG","TC","Height")){
    fill <- "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_"
  }else{
    fill <- "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_binary_"
  }
  dat <- read.csv(paste0(fill,trait,".regenie"), sep="")
  colnames(dat) <- c("CHR","BP","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
  dat$P <- 10^(-1*dat$LOG10P)
  
  Z <- as.numeric(dat$BETA/dat$SE)
  P <- as.numeric(dat$P)
  A1 <- dat$REF
  A2 <- dat$ALT
  N <- as.numeric(dat$N)
  
  SNP_ID <- dat$ID
  sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N)
  fwrite(sumstats,file = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDSC/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
  
  system(paste0("munge_sumstats.py",
                " --sumstats /data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDSC/sumstats_Overall",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDSC/sum_align",
                " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
                " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
                " --maf-min 0.01"))
  
  system(paste0("ldsc.py",
                " --h2 /data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDSC/sum_align.sumstats.gz",
                " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
                " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDSC/",trait,"_LDSC_Imputed_Results"))
}