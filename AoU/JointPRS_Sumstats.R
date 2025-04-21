rm(list = ls())

snpinfo_mult_1kg_hm3 <- read.delim("/data/williamsjacr/PRSCSx_LD/snpinfo_mult_1kg_hm3")

library(stringr)
library(data.table)
library(dplyr)

SNP_GRCh37_38_match_update <- readRDS("/data/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")
SNP_GRCh37_38_match_update$unique_id1 <- paste0("chr",SNP_GRCh37_38_match_update$chr,":",SNP_GRCh37_38_match_update$pos38,":",
                                                SNP_GRCh37_38_match_update$allele1_38,":",SNP_GRCh37_38_match_update$allele2_38)
SNP_GRCh37_38_match_update$unique_id2 <- paste0("chr",SNP_GRCh37_38_match_update$chr,":",SNP_GRCh37_38_match_update$pos38,":",
                                                SNP_GRCh37_38_match_update$allele2_38,":",SNP_GRCh37_38_match_update$allele1_38)

SNP_GRCh37_38_match_update$unique_id1 <- toupper(SNP_GRCh37_38_match_update$unique_id1)
SNP_GRCh37_38_match_update$unique_id1 <- gsub(" ","",SNP_GRCh37_38_match_update$unique_id1)
SNP_GRCh37_38_match_update$unique_id2 <- toupper(SNP_GRCh37_38_match_update$unique_id2)
SNP_GRCh37_38_match_update$unique_id2 <- gsub(" ","",SNP_GRCh37_38_match_update$unique_id2)

for(anc in c("EUR","AFR","AMR")){
  for(trait in c("BMI","LDL","HDL","Height","logTG","TC")){
    data <- fread(paste0("/data/williamsjacr/Cleaned_AoU_SumStats/",anc,"_",trait,"_GWAS_SumStats_Cleaned"))
    data <- as.data.frame(data)
    data <- data[,c("CHR","SNP","BP","REF","ALT","OBS_CT","BETA","SE","P","SNP")]
    colnames(data) <- c("CHR","SNP","BP","REF","ALT","OBS_CT","BETA","SE","P","rs_ID")
    data$SNP <- toupper(data$SNP)
    data$SNP <- gsub(" ","",data$SNP)
    
    data <- left_join(data,SNP_GRCh37_38_match_update[,c("unique_id1","rsid")],by = c("SNP"="unique_id1"))
    data$rs_ID[!is.na(data$rsid)] <- data$rsid[!is.na(data$rsid)]
    data <- subset(data,select = -c(rsid))
    data <- left_join(data,SNP_GRCh37_38_match_update[,c("unique_id2","rsid")],by = c("SNP"="unique_id2"))
    data$rs_ID[!is.na(data$rsid)] <- data$rsid[!is.na(data$rsid)]
    data <- data[,c("CHR","SNP","BP","REF","ALT","OBS_CT","BETA","SE","P","rs_ID")]
    data$SNP <- data$rs_ID
    doubles <- names(which(table(data$rs_ID) > 1))
    data <- data[!(data$SNP %in% doubles),]

    median_samplesize <- median(data$OBS_CT)
    data <- data[str_detect(data$SNP,"rs"),]
    data <- data[data$SNP %in% snpinfo_mult_1kg_hm3$SNP,]
    data <- data[!is.na(data$BETA),]
    JointPRS_SumStats <- data.frame(SNP = data$SNP,A1 = data$ALT,A2 = data$REF,BETA = data$BETA,P = data$P)
    write.table(JointPRS_SumStats,file = paste0("/data/williamsjacr/AoU_JointPRS/GWASSumStats/",anc,"_",trait,"_GWAS_SumStats_Cleaned.txt"),row.names = F,col.names = T,quote=F)
    BIM <- data.frame(CHR = data$CHR,SNP = data$SNP,v3 = 0,BP = data$BP,A1 = data$ALT,A2 = data$REF)
    write.table(BIM,file = paste0("/data/williamsjacr/AoU_JointPRS/GWASSumStats/",anc,"_",trait,".bim"),row.names = F,col.names = F,quote=F)
    save(median_samplesize,file = paste0("/data/williamsjacr/AoU_JointPRS/GWASSumStats/",anc,"_",trait,"_samplesize.RData"))
  }
}