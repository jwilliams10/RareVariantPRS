# %%writefile score_task.R
# 
# tasks <- data.frame(check.names = FALSE)
# 
# for(anc in c("AFR","AMR","EUR")){
#   tasks <- rbind(tasks, data.frame(
#     '--input BED_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.bed",
#     '--input BIM_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.bim",
#     '--input FAM_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.fam",
#     '--input all_phenotypes_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/all_phenotypes.csv",
#     '--input all_train_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/All_Train.txt",
#     '--input R_Script'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Scripts/reference_data.R",
#     '--env anc'=anc,
#     '--output-recursive OUTPUT_PATH'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants",
#     check.names = FALSE
#   ))
# }
# 
# 
# write.table(tasks, 
#             file="score_task.txt", 
#             row.names=F, col.names=T, 
#             sep='\t', quote=F)
# 
# %%writefile reference_data.R
rm(list = ls())

BED_file <- commandArgs(TRUE)[1]
BED_file <- gsub(".bed","",BED_file)
print(BED_file)

all_phenotypes_file <- commandArgs(TRUE)[2]
print(all_phenotypes_file)

all_train_file <- commandArgs(TRUE)[3]
print(all_phenotypes_file)

anc <- commandArgs(TRUE)[4]
print(anc)

OUTPUT_PATH <- commandArgs(TRUE)[5]
print(OUTPUT_PATH)

library(data.table)

bed_ids <- fread(paste0(BED_file,".fam"))
bed_ids <- as.data.frame(bed_ids)
bed_ids <- bed_ids[,1:2]

all_phenotypes <- read.csv(all_phenotypes_file)

tmp <- read.delim(all_train_file)

train_ids <- bed_ids[bed_ids[,2] %in% tmp$IID,]

tmp <- tmp[tmp$IID %in% all_phenotypes$IID[all_phenotypes$ancestry == anc],]

train_ids_anc <- train_ids[train_ids[,2] %in% tmp$IID,]

train_ids_ref <- train_ids_anc[sample(1:nrow(train_ids_anc),3000,replace = FALSE),]

write.table(train_ids_ref,"reference.txt",row.names = FALSE,col.names = FALSE)

system(paste0("plink2 --bfile ",BED_file," --keep reference.txt --make-bed --out ",OUTPUT_PATH,"/all_chr_reference_",anc))
