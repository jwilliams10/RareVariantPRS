rm(list = ls())

ukb_hm3_mega_bim <- read.delim("/data/BB_Bioinformatics/ProjectData/UKB/Genotypes/ukb_hm3_mega.bim", header=FALSE)

snp_list_chr_22 <- ukb_hm3_mega_bim[ukb_hm3_mega_bim[,1] == 22,2]

write.table(snp_list_chr_22,"/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/snp_list.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

sampleids <- read.table("/data/williamsjacr/UKB_WES_lipids/Data/sampleids.txt", quote="\"", comment.char="")

unrels_nRandomSNPs_0 <- read.table("/data/BB_Bioinformatics/ProjectData/UKB/phenotypes/unrels_nRandomSNPs_0.unrels", quote="\"", comment.char="")

load("/data/BB_Bioinformatics/ProjectData/UKB_WES_lipids/tmp.LDL.20211014.Rdata")
phenotype <- tmp.LDL
rm(tmp.LDL)

sampleids <- sampleids[sampleids[,1] %in% phenotype$userId,,drop = FALSE]
sampleids <- sampleids[sampleids[,1] %in% unrels_nRandomSNPs_0$V1,,drop = FALSE]

phenotype$FID <- 0

phenotype <- phenotype[,c("userId","FID","LDLadj.norm","age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]

colnames(phenotype) <- c("IID",colnames(phenotype)[-1])
phenotype[,c(4:5,7:16)] <- scale(phenotype[,c(4:5,7:16)]) 

phenotype <- phenotype[phenotype$IID %in% sampleids[,1],]

save(phenotype,file = "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sample_phenotype.RData")

write.table(cbind(sampleids,sampleids),"/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_common.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sampleids,"/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_rare.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)


#### add mean imputations
system("/data/williamsjacr/software/plink2 --bfile /data/BB_Bioinformatics/ProjectData/UKB/Genotypes/ukb_hm3_mega --extract /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/snp_list.txt --geno 0.1 --keep /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_common.txt --make-bed --out /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common")

