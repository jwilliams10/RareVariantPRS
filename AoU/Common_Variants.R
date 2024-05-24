# %%writefile score_task.R
# 
# tasks <- data.frame(check.names = FALSE)
# 
# for (chrom in 1:22) {
#   tasks <- rbind(tasks, data.frame(
#     '--env CHR'=chrom,
#     '--input BED_File'=paste0("gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr",chrom,".bed"),
#     '--input BIM_File'=paste0("gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr",chrom,".bim"),
#     '--input FAM_File'=paste0("gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/plink_bed/acaf_threshold.chr",chrom,".fam"),
#     '--input Ancestry_File'="gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv",
#     '--input R_Script'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Scripts/Common_Variants.R",
#     '--output-recursive OUTPUT_PATH'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants",
#     # '--output out_file'=paste0("gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/test",chrom,".Rdata"),
#     check.names = FALSE
#   ))
# }
# 
# write.table(tasks, 
#             file="score_task.txt", 
#             row.names=F, col.names=T, 
#             sep='\t', quote=F)


rm(list = ls())

chr <- as.numeric(commandArgs(TRUE)[1])
print(chr)

BED_file <- commandArgs(TRUE)[2]
BED_file <- gsub(".bed","",BED_file)
print(BED_file)

Ancestry_File <- commandArgs(TRUE)[5]
print(Ancestry_File)

OUTPUT_PATH <- commandArgs(TRUE)[6]
print(OUTPUT_PATH)


ancestries <- read.delim(Ancestry_File)
ancestries <- data.frame(IID = as.numeric(ancestries[,c("research_id")]),ancestry = toupper(ancestries[,c("ancestry_pred")]))


AFR_Keep <- data.frame(FID = 0,IID = ancestries$IID[ancestries$ancestry == "AFR"])
AMR_Keep <- data.frame(FID = 0,IID = ancestries$IID[ancestries$ancestry == "AMR"])
EUR_Keep <- data.frame(FID = 0,IID = ancestries$IID[ancestries$ancestry == "EUR"])


write.table(AFR_Keep,"AFR_Keep.txt",row.names = FALSE,col.names = FALSE)
write.table(AMR_Keep,"AMR_Keep.txt",row.names = FALSE,col.names = FALSE)
write.table(EUR_Keep,"EUR_Keep.txt",row.names = FALSE,col.names = FALSE)

system(paste0("plink --bfile ",BED_file," --keep AFR_Keep.txt --freq --out chr",chr,"_freq_AFR"),intern = TRUE)
system(paste0("plink --bfile ",BED_file," --keep AMR_Keep.txt --freq --out chr",chr,"_freq_AMR"),intern = TRUE)
system(paste0("plink --bfile ",BED_file," --keep EUR_Keep.txt --freq --out chr",chr,"_freq_EUR"),intern = TRUE)

a <- read.csv(paste0("chr",chr,"_freq_AFR.frq"),sep = "")
b <- read.csv(paste0("chr",chr,"_freq_AMR.frq"),sep = "")
c <- read.csv(paste0("chr",chr,"_freq_EUR.frq"),sep = "")

exclude_list <- a$SNP[is.na(a$MAF) | is.na(b$MAF) | is.na(c$MAF) | a$MAF <= 0.01 | b$MAF <= 0.01 | c$MAF <= 0.01]
write.table(exclude_list,"exclude_list.txt",row.names = F,col.names = F,quote=F)

system(paste0("plink --bfile ",BED_file," --keep AFR_Keep.txt --exclude exclude_list.txt --make-bed --out ",OUTPUT_PATH,"/AFR_chr",chr,"_corrected"),intern = TRUE)
system(paste0("plink --bfile ",BED_file," --keep AMR_Keep.txt --exclude exclude_list.txt --make-bed --out ",OUTPUT_PATH,"/AMR_chr",chr,"_corrected"),intern = TRUE)
system(paste0("plink --bfile ",BED_file," --keep EUR_Keep.txt --exclude exclude_list.txt --make-bed --out ",OUTPUT_PATH,"/EUR_chr",chr,"_corrected"),intern = TRUE)