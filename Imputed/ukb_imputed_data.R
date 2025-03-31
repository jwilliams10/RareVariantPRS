rm(list = ls())
library(SeqArray)

All_Train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
All_Tune <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
All_Validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")

phenotypes <- rbind(All_Train,All_Tune,All_Validation)

fam_file <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/Genotypes/ukb_hm3_mega.fam", quote="\"", comment.char="")

gds_file <- seqOpen("/data/williamsjacr/UKB_WES_Full_Processed_Data/gds/chr21.gds") 
gds_ids <- seqGetData(gds_file,"sample.id")
seqClose(gds_file)

intersection_ids <- intersect(phenotypes$IID,fam_file$V1)
intersection_ids <- intersect(intersection_ids,gds_ids)
intersection_ids <- as.numeric(intersection_ids)

write.table(cbind(intersection_ids,intersection_ids),"/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/intersection_ids.txt",row.names = FALSE,col.names = FALSE)

system("/data/williamsjacr/software/plink2 --bfile /data/BB_Bioinformatics/ProjectData/UKB/Genotypes/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/intersection_ids.txt --maf 0.01 --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega")

fam_file <- read.table("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega.fam", quote="\"", comment.char="")
fam_file$V1 <- 0
write.table(fam_file,file = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega.fam",row.names = FALSE,col.names = FALSE)

all_chr <- read.table("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega.bim", header=FALSE)
samps <- sample(1:nrow(all_chr),500000,replace = FALSE)
samps <- samps[order(samps)]
all_chr <- all_chr[samps,2]
write.table(all_chr,file = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/regenie_extract_snps.txt",row.names = F,col.names = F,quote=F)

system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --extract /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/regenie_extract_snps.txt --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1")