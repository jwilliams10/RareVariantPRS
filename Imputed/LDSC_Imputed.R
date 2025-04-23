rm(list = ls())

library(data.table)
library(dplyr)

for(trait in c("BMI","LDL_statin_adj","HDL","logTG","TC_statin_adj","Height")){
  dat <- read.delim(paste0("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/GWAS_SummaryStatistics/",trait,"_sumstats.",trait,".glm.linear"), header=FALSE, comment.char="#")
  colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
  dat <- dat[dat$TEST == "ADD",]
  
  Z <- as.numeric(dat$BETA/dat$SE)
  P <- as.numeric(dat$P)
  A1 <- dat$REF
  A2 <- dat$ALT
  N <- as.numeric(dat$OBS_CT)
  
  SNP_ID <- dat$ID
  sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N)
  fwrite(sumstats,file = "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
  
  system(paste0("munge_sumstats.py",
                " --sumstats /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sumstats_Overall",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sum_align",
                " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
                " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
                " --maf-min 0.01"))
  
  system(paste0("ldsc.py",
                " --h2 /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sum_align.sumstats.gz",
                " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
                " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Imputed_Results"))
}