rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/FAVORdatabase_chrsplit.csv"

### xsv directory
xsv <- "~/.cargo/bin/xsv"

### output (Step 2 and Step 3)
output_path <- "/data/williamsjacr/UKB_WES_lipids/Data/agds/"

### DB file
DB_path <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/"

### anno channel (subset)
anno_colnum <- c(1,8,9:12,14,16,19,23,25:36)

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### chromosome number
## annotate (seperate)
DB_info <- read.csv(file_DBsplit,header=TRUE)
chr_splitnum <- sum(DB_info$Chr==chr)

for(kk in 1:chr_splitnum)
{
  print(kk)
  
  system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/VarInfo_train_chr",chr,"_",kk,".csv variant_vcf ",DB_path,"/chr",chr,"_",kk,".csv > ",output_path,"chr",chr,"/Anno_chr_train",chr,"_",kk,".csv"))
}

## merge info
Anno <- paste0(output_path,"chr",chr,"/Anno_chr_train",chr,"_",seq(1:chr_splitnum),".csv ")
merge_command <- paste0(xsv," cat rows ",Anno[1])

for(kk in 2:chr_splitnum)
{
  merge_command <- paste0(merge_command,Anno[kk])
}

merge_command <- paste0(merge_command,"> ",output_path,"chr",chr,"/Anno_chr_train",chr,".csv")

system(merge_command)

## subset
anno_colnum_xsv <- c()
for(kk in 1:(length(anno_colnum)-1))
{
  anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[kk],",")
}
anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[length(anno_colnum)])

system(paste0(xsv," select ",anno_colnum_xsv," ",output_path,"chr",chr,"/Anno_chr_train",chr,".csv > ",output_path,"chr",chr,"/Anno_chr_train",chr,"_STAARpipeline.csv"))














rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/FAVORdatabase_chrsplit.csv"

### xsv directory
xsv <- "~/.cargo/bin/xsv"

### output (Step 2 and Step 3)
output_path <- "/data/williamsjacr/UKB_WES_lipids/Data/agds/"

### DB file
DB_path <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/"

### anno channel (subset)
anno_colnum <- c(1,8,9:12,14,16,19,23,25:36)

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### chromosome number
## annotate (seperate)
DB_info <- read.csv(file_DBsplit,header=TRUE)
chr_splitnum <- sum(DB_info$Chr==chr)

for(kk in 1:chr_splitnum)
{
  print(kk)
  
  system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/VarInfo_tune_chr",chr,"_",kk,".csv variant_vcf ",DB_path,"/chr",chr,"_",kk,".csv > ",output_path,"chr",chr,"/Anno_chr_tune",chr,"_",kk,".csv"))
}

## merge info
Anno <- paste0(output_path,"chr",chr,"/Anno_chr_tune",chr,"_",seq(1:chr_splitnum),".csv ")
merge_command <- paste0(xsv," cat rows ",Anno[1])

for(kk in 2:chr_splitnum)
{
  merge_command <- paste0(merge_command,Anno[kk])
}

merge_command <- paste0(merge_command,"> ",output_path,"chr",chr,"/Anno_chr_tune",chr,".csv")

system(merge_command)

## subset
anno_colnum_xsv <- c()
for(kk in 1:(length(anno_colnum)-1))
{
  anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[kk],",")
}
anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[length(anno_colnum)])

system(paste0(xsv," select ",anno_colnum_xsv," ",output_path,"chr",chr,"/Anno_chr_tune",chr,".csv > ",output_path,"chr",chr,"/Anno_chr_tune",chr,"_STAARpipeline.csv"))
















rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/FAVORdatabase_chrsplit.csv"

### xsv directory
xsv <- "~/.cargo/bin/xsv"

### output (Step 2 and Step 3)
output_path <- "/data/williamsjacr/UKB_WES_lipids/Data/agds/"

### DB file
DB_path <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/"

### anno channel (subset)
anno_colnum <- c(1,8,9:12,14,16,19,23,25:36)

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### chromosome number
## annotate (seperate)
DB_info <- read.csv(file_DBsplit,header=TRUE)
chr_splitnum <- sum(DB_info$Chr==chr)

for(kk in 1:chr_splitnum)
{
  print(kk)
  
  system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/VarInfo_validation_chr",chr,"_",kk,".csv variant_vcf ",DB_path,"/chr",chr,"_",kk,".csv > ",output_path,"chr",chr,"/Anno_chr_validation",chr,"_",kk,".csv"))
}

## merge info
Anno <- paste0(output_path,"chr",chr,"/Anno_chr_validation",chr,"_",seq(1:chr_splitnum),".csv ")
merge_command <- paste0(xsv," cat rows ",Anno[1])

for(kk in 2:chr_splitnum)
{
  merge_command <- paste0(merge_command,Anno[kk])
}

merge_command <- paste0(merge_command,"> ",output_path,"chr",chr,"/Anno_chr_validation",chr,".csv")

system(merge_command)

## subset
anno_colnum_xsv <- c()
for(kk in 1:(length(anno_colnum)-1))
{
  anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[kk],",")
}
anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[length(anno_colnum)])

system(paste0(xsv," select ",anno_colnum_xsv," ",output_path,"chr",chr,"/Anno_chr_validation",chr,".csv > ",output_path,"chr",chr,"/Anno_chr_validation",chr,"_STAARpipeline.csv"))




