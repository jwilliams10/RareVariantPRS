rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/FAVORdatabase_chrsplit.csv"
### Targeted GDS
dir_geno <- "/data/williamsjacr/UKB_WES_lipids/Data/gds/"
gds_file_name_1 <- "train_chr"
gds_file_name_2 <- ".gds"
### output
output_path <- "/data/williamsjacr/UKB_WES_lipids/Data/agds/"

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### make directory
system(paste0("mkdir ",output_path,"chr",chr))

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)

seqClose(genofile)

## Generate VarInfo
for(kk in 1:dim(DB_info)[1])
{
  print(kk)
  
  VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
  VarInfo <- data.frame(VarInfo)
  
  write.csv(VarInfo,paste0(output_path,"chr",chr,"/VarInfo_train_chr",chr,"_",kk,".csv"),quote=FALSE,row.names = FALSE)
}







rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/FAVORdatabase_chrsplit.csv"
### Targeted GDS
dir_geno <- "/data/williamsjacr/UKB_WES_lipids/Data/gds/"
gds_file_name_1 <- "tune_chr"
gds_file_name_2 <- ".gds"
### output
output_path <- "/data/williamsjacr/UKB_WES_lipids/Data/agds/"

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### make directory
system(paste0("mkdir ",output_path,"chr",chr))

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)

seqClose(genofile)

## Generate VarInfo
for(kk in 1:dim(DB_info)[1])
{
  print(kk)
  
  VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
  VarInfo <- data.frame(VarInfo)
  
  write.csv(VarInfo,paste0(output_path,"chr",chr,"/VarInfo_tune_chr",chr,"_",kk,".csv"),quote=FALSE,row.names = FALSE)
}







rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
file_DBsplit <- "/data/BB_Bioinformatics/FAVOR_Annotations/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB/FAVORdatabase_chrsplit.csv"
### Targeted GDS
dir_geno <- "/data/williamsjacr/UKB_WES_lipids/Data/gds/"
gds_file_name_1 <- "validation_chr"
gds_file_name_2 <- ".gds"
### output
output_path <- "/data/williamsjacr/UKB_WES_lipids/Data/agds/"

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function 
###########################################################################

### make directory
system(paste0("mkdir ",output_path,"chr",chr))

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)

seqClose(genofile)

## Generate VarInfo
for(kk in 1:dim(DB_info)[1])
{
  print(kk)
  
  VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
  VarInfo <- data.frame(VarInfo)
  
  write.csv(VarInfo,paste0(output_path,"chr",chr,"/VarInfo_validation_chr",chr,"_",kk,".csv"),quote=FALSE,row.names = FALSE)
}

