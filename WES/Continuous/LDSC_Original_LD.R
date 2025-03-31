rm(list = ls())

library(data.table)
library(dplyr)

SNPinfo <- readRDS("/data/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")

SNPinfo_38 <- data.frame(rsid = SNPinfo$rsid,unique_id = paste0(SNPinfo$chr,"_",SNPinfo$pos38,"_",SNPinfo$allele1_38,"_",SNPinfo$allele2_38))

rm(SNPinfo)

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  dat <- read.delim(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/GWAS_Summary_Statistics/",trait,"_sumstats.",trait,".glm.linear"), header=FALSE, comment.char="#")
  colnames(dat) <- c("CHR","BP","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
  dat <- dat[dat$TEST == "ADD",]
  
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
  N <- as.numeric(dat$OBS_CT)
  
  SNP_ID <- dat$ID
  sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N)
  fwrite(sumstats,file = "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
  
  system(paste0("munge_sumstats.py",
                " --sumstats /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sumstats_Overall",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sum_align",
                " --merge-alleles /data/BB_Bioinformatics/ProjectData/LDSC_reference/original_EUR_HM3/eur_w_ld_chr/w_hm3.snplist",
                " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
                " --maf-min 0.01"))
  
  system(paste0("ldsc.py",
                " --h2 /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/sum_align.sumstats.gz",
                " --ref-ld-chr /data/BB_Bioinformatics/ProjectData/LDSC_reference/original_EUR_HM3/eur_w_ld_chr/",
                " --w-ld-chr /data/BB_Bioinformatics/ProjectData/LDSC_reference/original_EUR_HM3/eur_w_ld_chr/",
                " --out /data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDSC/",trait,"_LDSC_Original_Results"))
}