rm(list = ls())

library(data.table)
library(dplyr)

SNPinfo <- readRDS("/data/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")

SNPinfo_38 <- data.frame(rsid = SNPinfo$rsid,unique_id = paste0(SNPinfo$chr,"_",SNPinfo$pos38,"_",SNPinfo$allele1_38,"_",SNPinfo$allele2_38))

rm(SNPinfo)

for(trait in c("Breast","Prostate","CAD","T2D","Asthma")){
  dat <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step2_binary_",trait,".regenie"), sep="")
  colnames(dat) <- c("CHR","BP","ID","REF","ALT","A1_FREQ","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
  dat$P <- 10^(-1*dat$LOG10P)
  
  dat$unique_id1 <- paste0(dat$CHR,"_",dat$BP,"_",
                           dat$REF,"_",dat$ALT)
  dat$unique_id2 <- paste0(dat$CHR,"_",dat$BP,"_",
                           dat$ALT,"_",dat$REF)
  
  dat <- left_join(dat,SNPinfo_38,by = c("unique_id1" = "unique_id"))
  dat$ID[!is.na(dat$rsid)] <- dat$rsid[!is.na(dat$rsid)]
  
  dat <- subset(dat,select = -c(rsid))
  
  dat <- left_join(dat,SNPinfo_38,by = c("unique_id2" = "unique_id"))
  dat$ID[!is.na(dat$rsid)] <- dat$rsid[!is.na(dat$rsid)]
  
  Z <- as.numeric(dat$BETA/dat$SE)
  P <- as.numeric(dat$P)
  A1 <- dat$REF
  A2 <- dat$ALT
  N <- as.numeric(dat$N)
  
  SNP_ID <- dat$ID
  sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N)
  fwrite(sumstats,file = "/data/williamsjacr/UKB_WES_Phenotypes/Binary/LDSC/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
  
  system(paste0("munge_sumstats.py",
                " --sumstats /data/williamsjacr/UKB_WES_Phenotypes/Binary/LDSC/sumstats_Overall",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/LDSC/sum_align",
                " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
                " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
                " --maf-min 0.01"))
  
  system(paste0("ldsc.py",
                " --h2 /data/williamsjacr/UKB_WES_Phenotypes/Binary/LDSC/sum_align.sumstats.gz",
                " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
                " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/LDSC/",trait,"_LDSC_Results"))
}