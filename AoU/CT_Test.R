# %%writefile score_task.R
# 
# tasks <- data.frame(check.names = FALSE)
# 
# anc <- "EUR"
# trait <- "HDL"
# 
# tasks <- rbind(tasks, data.frame(
#   '--input BED_Full_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.bed",
#   '--input BIM_Full_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.bim",
#   '--input FAM_Full_File'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr.fam",
#   '--input BED_Ref_File'=paste0("gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr_reference_",anc,".bed"),
#   '--input BIM_Ref_File'=paste0("gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr_reference_",anc,".bim"),
#   '--input FAM_Ref_File'=paste0("gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants/all_chr_reference_",anc,".fam"),
#   '--input sumstats'=paste0("gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/Continuous/GWAS_SummaryStats/",trait,"_sumstats_",anc,".",trait,".glm.linear"),
#   '--input all_phenotypes_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/all_phenotypes.csv",
#   '--input all_train_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/All_Train.txt",
#   '--input all_tune_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/All_Tune.txt",
#   '--input all_valid_file'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/All_Validation.txt",
#   '--input R_Script'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Scripts/CT_Test.R",
#   '--env anc'=anc,
#   '--env trait'=trait,
#   '--output-recursive OUTPUT_PATH'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/Continuous/Test_CT",
#   check.names = FALSE
# ))
# 
# 
# 
# write.table(tasks, 
#             file="score_task.txt", 
#             row.names=F, col.names=T, 
#             sep='\t', quote=F)
# 
# %%writefile CT_Test.R
rm(list = ls())

BED_Full_File <- commandArgs(TRUE)[1]
BED_Full_File <- gsub(".bed","",BED_Full_File)
print(BED_Full_File)

BED_Ref_File <- commandArgs(TRUE)[2]
BED_Ref_File <- gsub(".bed","",BED_Ref_File)
print(BED_Ref_File)

sumstats <- commandArgs(TRUE)[3]
print(sumstats)

all_phenotypes_file <- commandArgs(TRUE)[4]
print(all_phenotypes_file)

all_train_file <- commandArgs(TRUE)[5]
print(all_train_file)

all_tune_file <- commandArgs(TRUE)[6]
print(all_tune_file)

all_valid_file <- commandArgs(TRUE)[7]
print(all_valid_file)

anc <- commandArgs(TRUE)[8]
print(anc)

trait <- commandArgs(TRUE)[9]
print(trait)

OUTPUT_PATH <- commandArgs(TRUE)[10]
print(OUTPUT_PATH)

library(data.table)
library(dplyr)

bed_ids <- fread(paste0(BED_Full_File,".fam"))
bed_ids <- as.data.frame(bed_ids)
bed_ids <- bed_ids[,1:2]

all_phenotypes <- read.csv(all_phenotypes_file)



tmp <- read.delim(all_train_file)

train_ids <- bed_ids[bed_ids[,2] %in% tmp$IID,]

tmp <- tmp[tmp$IID %in% all_phenotypes$IID[all_phenotypes$ancestry == anc],]

train_ids_anc <- train_ids[train_ids[,2] %in% tmp$IID,]
write.table(train_ids_anc,"train.txt",row.names = FALSE,col.names = FALSE)



tmp <- read.delim(all_tune_file)

tune_ids <- bed_ids[bed_ids[,2] %in% tmp$IID,]

tmp <- tmp[tmp$IID %in% all_phenotypes$IID[all_phenotypes$ancestry == anc],]

tune_ids_anc <- tune_ids[tune_ids[,2] %in% tmp$IID,]
write.table(tune_ids_anc,"tune.txt",row.names = FALSE,col.names = FALSE)



tmp <- read.delim(all_valid_file)

valid_ids <- bed_ids[bed_ids[,2] %in% tmp$IID,]

tmp <- tmp[tmp$IID %in% all_phenotypes$IID[all_phenotypes$ancestry == anc],]

valid_ids_anc <- valid_ids[valid_ids[,2] %in% tmp$IID,]
write.table(valid_ids_anc,"validation.txt",row.names = FALSE,col.names = FALSE)







dat <- read.delim(sumstats, header=FALSE, comment.char="#")
colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
dat <- dat[dat$TEST == "ADD",]

dat <- dat[,c("CHROM","ID","REF","POS","A1","BETA","P")]
colnames(dat) <- c("CHR","SNP","REF","BP","A1","BETA","P")

write.table(dat,file = paste0(trait,"_assoc.txt"),col.names = T,row.names = F,quote=F)

pthr <- 1
r2thr <- 0.1
kbpthr <- 500

# Run plink 1.9 with summary statistics data with --bfile being the 1000 genomes reference file, --clump being the summary statistics file, and this is written to temp.dir/LD_clump
system(paste0("plink --bfile ",BED_Ref_File," --clump ",trait,"_assoc.txt --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",trait,"_LDclump"))

################################################

## This is the beginning of the thresholding step of clumping and threshold

#p-value thresholds
LD <- as.data.frame(fread(paste0(trait,"_LDclump.clumped")))

# grab the index SNP for each clump
clump.snp <- LD[,3,drop=F]
# join against the summary data, this is now a dataset with n = number of index SNPs
prs.all <- left_join(clump.snp,dat)

## Write Coefficients of index SNPs to use later
prs.file <- prs.all[,c("SNP","A1","BETA")]
write.table(prs.file,file = paste0(trait,"_prs_coeff"),col.names = T,row.names = F,quote=F)

# Write p-values to file to use later
p.value.file <- prs.all[,c("SNP","P")]
write.table(p.value.file,file = paste0(trait,"_p_value"),col.names = T,row.names = F,quote=F)

beta_final <- inner_join(prs.file,p.value.file)

pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
n_pthres <- length(pthres)
## Write a file with the p-value thresholds 
q_range <- data.frame(a = paste0("p_value_",1:length(pthres)),b = 0,c = pthres)

write.table(q_range,file = "q_range_file",row.names = F,col.names = F,quote=F)

#PRS = G*beta/(2*number of SNPs) #column header is SCORE_AVG
#PRS = G*beta
# Literally just multiply the score file (weighted or unweighted coefficients) by the G matrix, q-score-range is only for C + T, for LD pred score file would be weight coefficients 

system(paste0("plink2 --q-score-range q_range_file ",trait,"_p_value header --threads 2 --score ",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile ",BED_Full_File," --keep train.txt --out ",trait,"_prs_train"))
system(paste0("plink2 --q-score-range q_range_file ",trait,"_p_value header --threads 2 --score ",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile ",BED_Full_File," --keep tune.txt --out ",trait,"_prs_tune"))
system(paste0("plink2 --q-score-range q_range_file ",trait,"_p_value header --threads 2 --score ",trait,"_prs_coeff cols=+scoresums,-scoreavgs header no-mean-imputation  --bfile ",BED_Full_File," --keep validation.txt --out ",trait,"_prs_validation"))
#########################################################################

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0(trait,"_prs_train.p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_train <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_train)[2] <- "IID"

write.table(prs_mat_train,file = paste0(OUTPUT_PATH,"/",trait,"_prs_all_train.txt"),row.names = F)

prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("",trait,"_prs_tune.p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_tune <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_tune)[2] <- "IID"

write.table(prs_mat_tune,file = paste0(OUTPUT_PATH,"/",trait,"_prs_all_tune.txt"),row.names = F)

### Merge all the outputted files from the previous command into one large data file
prs_list <- list()
temp <- 1
for(k in 1:length(pthres)){
  
  prs_temp <- fread(paste0("",trait,"_prs_validation.p_value_",k,".sscore"))
  # times (2*number of SNPs)
  prs_list[[temp]] <- prs_temp[,5:ncol(prs_temp)]
  
  colnames(prs_list[[temp]]) <- paste0("p_value_",k)
  temp <- temp + 1
  
}
prs_mat_validation <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
colnames(prs_mat_validation)[2] <- "IID"

write.table(prs_mat_validation,file = paste0(OUTPUT_PATH,"/",trait,"_prs_all_validation.txt"),row.names = F)