rm(list = ls())

ukb_hm3_mega_bim <- read.delim("/data/BB_Bioinformatics/ProjectData/UKB/Genotypes/ukb_hm3_mega.bim", header=FALSE)

snp_list_chr_22 <- ukb_hm3_mega_bim[ukb_hm3_mega_bim[,1] == 22,2]

write.table(snp_list_chr_22,"/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/snp_list.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

sampleids <- read.table("/data/williamsjacr/UKB_WES_Full_Processed_Data/sampleids.txt", quote="\"", comment.char="")

unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

sampleids <- sampleids[sampleids[,1] %in% unrels_nRandomSNPs_0$V1,,drop = FALSE]
sampleids <- sampleids[sampleids[,1] %in% ukb_pheno$IID,,drop = FALSE]

phenotype <- data.frame(IID = sampleids,FID = 0,Y = rnorm(nrow(sampleids)))

save(phenotype,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sample_phenotype.RData")

write.table(cbind(sampleids,sampleids),"/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_common.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sampleids,"/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_rare.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)


#### add mean imputations
system("/data/williamsjacr/software/plink2 --bfile /data/BB_Bioinformatics/ProjectData/UKB/Genotypes/ukb_hm3_mega --extract /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/snp_list.txt --geno 0.1 --keep /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_common.txt --make-bed --out /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common")

chr22_filtered_common <- read.delim("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam", header=FALSE)
chr22_filtered_common$V1 <- 0
write.table(chr22_filtered_common,col.names = FALSE,row.names = FALSE,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam")
