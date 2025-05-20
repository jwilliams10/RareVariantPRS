rm(list = ls())

library(data.table)
library(dplyr)

SNPinfo <- readRDS("/data/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")

SNPinfo_38 <- data.frame(rsid = SNPinfo$rsid,unique_id = paste0(SNPinfo$chr,"_",SNPinfo$pos38,"_",SNPinfo$allele1_38,"_",SNPinfo$allele2_38))

rm(SNPinfo)

for(anc in c("EUR","AMR","AFR")){
  for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
    
    dat <- fread(paste0("/data/williamsjacr/Cleaned_AoU_SumStats/",anc,"_",trait,"_GWAS_SumStats_Cleaned"))
    dat <- as.data.frame(dat)
    dat <- dat[,c("CHR","SNP","BP","REF","ALT","OBS_CT","BETA","SE","P","SNP")]
    colnames(dat) <- c("CHR","SNP","BP","REF","ALT","OBS_CT","BETA","SE","P","rs_ID")
    dat$SNP <- toupper(dat$SNP)
    dat$SNP <- gsub(" ","",dat$SNP)
    
    dat$unique_id1 <- paste0(dat$CHR,"_",dat$BP,"_",
                             dat$REF,"_",dat$ALT)
    dat$unique_id2 <- paste0(dat$CHR,"_",dat$BP,"_",
                             dat$ALT,"_",dat$REF)
    
    dat <- left_join(dat,SNPinfo_38,by = c("unique_id1" = "unique_id"))
    dat$rs_ID[!is.na(dat$rsid)] <- dat$rsid[!is.na(dat$rsid)]
    
    dat <- subset(dat,select = -c(rsid))
    
    dat <- left_join(dat,SNPinfo_38,by = c("unique_id2" = "unique_id"))
    dat$rs_ID[!is.na(dat$rsid)] <- dat$rsid[!is.na(dat$rsid)]
    
    Z <- as.numeric(dat$BETA/dat$SE)
    P <- as.numeric(dat$P)
    A1 <- dat$REF
    A2 <- dat$ALT
    N <- as.numeric(dat$OBS_CT)
    
    SNP_ID <- dat$rs_ID
    sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N)
    sumstats <- sumstats[sumstats$Z^2 < 80,]
    fwrite(sumstats,file = "/data/williamsjacr/Cleaned_AoU_SumStats/LDSC/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
    
    system(paste0("munge_sumstats.py",
                  " --sumstats /data/williamsjacr/Cleaned_AoU_SumStats/LDSC/sumstats_Overall",
                  " --out /data/williamsjacr/Cleaned_AoU_SumStats/LDSC/sum_align",
                  " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/",anc,"/",anc,"_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
                  " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
                  " --maf-min 0.01"))
    
    system(paste0("ldsc.py",
                  " --h2 /data/williamsjacr/Cleaned_AoU_SumStats/LDSC/sum_align.sumstats.gz",
                  " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/",anc,"/",anc,"_hm3_mega_ldsc/",
                  " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/",anc,"/",anc,"_hm3_mega_ldsc/",
                  " --out /data/williamsjacr/Cleaned_AoU_SumStats/LDSC/",anc,"_",trait,"_LDSC_Results"))
  } 
}