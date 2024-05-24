# %%writefile score_task.R
# 
# tasks <- data.frame(check.names = FALSE)
# 
# for(anc in c("AFR","AMR","EUR")){
#   for(trait in c("BMI","LDL","HDL","TC","logTG","Height")){
#     tasks <- rbind(tasks, data.frame(
#       '--input BED_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.bed",
#       '--input BIM_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.bim",
#       '--input FAM_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.fam",
#       '--input all_phenotypes_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/all_phenotypes.csv",
#       '--input all_train_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/All_Train.txt",
#       '--input R_Script'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Scripts/GWAS_SumStats.R",
#       '--env anc'=anc,
#       '--env trait'=trait,
#       '--output-recursive OUTPUT_PATH'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/Continuous/GWAS_SummaryStats",
#       check.names = FALSE
#     ))
#   }
# }
# 
# 
# write.table(tasks, 
#             file="score_task.txt", 
#             row.names=F, col.names=T, 
#             sep='\t', quote=F)
# 
# 
# %%writefile GWAS_SumStats.R
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

trait <- commandArgs(TRUE)[5]
print(trait)

OUTPUT_PATH <- commandArgs(TRUE)[6]
print(OUTPUT_PATH)

all_phenotypes <- read.csv(all_phenotypes_file)

tmp <- read.delim(all_train_file)
write.table(tmp[tmp$IID %in% all_phenotypes$IID[all_phenotypes$ancestry == anc],],file = paste0("All_Train_",anc,".txt"),sep = '\t',row.names = FALSE,quote = FALSE)

system(paste0("plink2 --bfile ",BED_file," --pheno All_Train_",anc,".txt --pheno-name ",trait," --linear --covar All_Train_",anc,".txt --covar-name age, age2, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 --vif 999 --out ",OUTPUT_PATH,"/",trait,"_sumstats_",anc))

