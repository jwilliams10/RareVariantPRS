rm(list = ls())

snpinfo_mult_1kg_hm3 <- read.delim("/data/williamsjacr/PRSCSx_LD/snpinfo_mult_1kg_hm3")

library(stringr)

for(anc in c("EUR","AFR","AMR")){
  for(trait in c("BMI","LDL","HDL","Height","logTG","TC")){
    data <- read.csv(paste0("/data/williamsjacr/Cleaned_AoU_SumStats/",anc,"_",trait,"_GWAS_SumStats_Cleaned.csv"))
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