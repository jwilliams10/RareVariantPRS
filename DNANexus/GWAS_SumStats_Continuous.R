rm(list = ls())

# for array in 1 2 3 4 5 6;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Clean_Data/all_chr.bed -iin=UKB_PRS:JW/Clean_Data/all_chr.bim -iin=UKB_PRS:JW/Clean_Data/all_chr.fam -iin=UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/GWAS_SumStats_Continuous.R -icmd="Rscript GWAS_SumStats_Continuous.R ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/ --instance-type mem1_ssd1_v2_x36
# done

a <- as.numeric(commandArgs(TRUE)[1])

trait <- c("BMI","Height","LDL","HDL","logTG","TC")[a]

All_Train <- read.delim("All_Train.txt")
All_Train <- All_Train[!is.na(All_Train[,trait]),c(2,1,3:ncol(All_Train))]
All_Train$FID <- 0
write.table(All_Train,file = "All_Train_REGENIE.txt",sep = '\t',row.names = FALSE,quote = FALSE)
print(nrow(All_Train))

set.seed(1330)

all_chr <- read.table("all_chr.bim", header=FALSE)
samps <- sample(1:nrow(all_chr),900000,replace = FALSE)
samps <- samps[order(samps)]
all_chr <- all_chr[samps,2]
write.table(all_chr,file = "extract_snps",row.names = F,col.names = F,quote=F)

system("plink2 --bfile all_chr --extract extract_snps --make-bed --out subsetted_allchr")

system("rm extract_snps")

subsetted_allchr <- read.table("subsetted_allchr.fam", quote="\"", comment.char="")
subsetted_allchr[,1] <- 0  
write.table(subsetted_allchr,"subsetted_allchr.fam",row.names = FALSE,col.names = FALSE)

allchr <- read.table("all_chr.fam", quote="\"", comment.char="")
allchr[,1] <- 0  
write.table(allchr,"all_chr.fam",row.names = FALSE,col.names = FALSE)

system(paste0("regenie --step 1 --bed subsetted_allchr --phenoFile All_Train_REGENIE.txt --phenoColList ",trait," --covarFile All_Train_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --qt --bsize 1000 --out regenie_step1_",trait))
system(paste0("regenie --step 2 --bed all_chr --phenoFile All_Train_REGENIE.txt --phenoColList ",trait," --covarFile All_Train_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred regenie_step1_",trait,"_pred.list --qt --bsize 400 --out regenie_step2_",trait))

system("rm all_chr.bed")
system("rm all_chr.fam")
system("rm all_chr.bim")

system("rm subsetted_allchr.bed")
system("rm subsetted_allchr.fam")
system("rm subsetted_allchr.bim")
system("rm All_Train.txt")
system("rm All_Train_REGENIE.txt")