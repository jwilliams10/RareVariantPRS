rm(list = ls())

dat <- NULL
for(i in 1:22){
  dat <- rbind(dat, read.delim(paste0("/data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/LDL_adj_norm__chr",i,"_train.LDLadj.norm.glm.linear"), header=FALSE, comment.char="#"))
}

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")

dat <- dat[,c("CHROM","ID","POS","A1","BETA","P")]
colnames(dat) <- c("CHR","SNP","BP","A1","BETA","P")

write.table(dat,file = "/data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/all_chr_assoc.txt",col.names = T,row.names = F,quote=F)

pthr <- 1
r2thr <- 0.1
kbpthr <- 500
# Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
system(paste0("/data/williamsjacr/software/plink --bfile ","/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr21_common_reference --clump /data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/all_chr_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ","/data/williamsjacr/UKB_WES_lipids/Data/Results/LDL/CT/LD_clump"))
