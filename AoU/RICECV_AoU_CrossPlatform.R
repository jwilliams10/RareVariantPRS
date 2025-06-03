rm(list = ls())

library(data.table)
library(dplyr)
library(boot)

trait <- c("Height","BMI","TC","HDL","LDL","logTG")[as.numeric(commandArgs(TRUE)[1])]

Score <- fread(paste0("/data/williamsjacr/AoU_Results/",trait,"_Final_Score"))
Score <- as.data.frame(Score)
name_split <- strsplit(Score$SNP,":")
CHR <- lapply(name_split,function(x){x[1]})
CHR <- as.numeric(gsub("chr","",CHR))
POS <- as.numeric(lapply(name_split,function(x){x[2]}))
A0 <- as.character(lapply(name_split,function(x){x[3]}))
A1 <- as.character(lapply(name_split,function(x){x[4]}))


SNPinfo <- readRDS("/data/BB_Bioinformatics/ProjectData/SNP_infor_GRCh37_38/SNP_GRCh37_38_match_update.rds")

SNPinfo_38 <- data.frame(rsid = SNPinfo$rsid,unique_id = paste0(SNPinfo$chr,"_",SNPinfo$pos38,"_",SNPinfo$allele1_38,"_",SNPinfo$allele2_38))

rm(SNPinfo)

Score$unique_id1 <- paste0(CHR,"_",POS,"_",A0,"_",A1)
Score$unique_id2 <- paste0(CHR,"_",POS,"_",A1,"_",A0)
Score$rs_ID <- NA

Score <- left_join(Score,SNPinfo_38,by = c("unique_id1" = "unique_id"))
Score$rs_ID[!is.na(Score$rsid)] <- Score$rsid[!is.na(Score$rsid)] 
Score <- subset(Score,select = -rsid)

Score <- left_join(Score,SNPinfo_38,by = c("unique_id2" = "unique_id"))
Score$rs_ID[!is.na(Score$rsid)] <- Score$rsid[!is.na(Score$rsid)] 
Score <- subset(Score,select = -rsid)

print(sum(is.na(Score$rs_ID))/nrow(Score))

Score$SNP <- Score$rs_ID
Score <- Score[,c("SNP","A1","BETA")]

write.table(Score,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Final_Score"),col.names = T,row.names = F,quote=F)

system(paste0("/data/williamsjacr/software/plink2 --score /data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Final_Score cols=+scoresums,-scoreavgs header no-mean-imputation --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --threads 1 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_RICECV_PRS"))