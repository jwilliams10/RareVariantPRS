rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### gds file
dir_geno <- "/data/williamsjacr/UKB_WES_Full_Processed_Data/gds/"
gds_file_name_1 <- "chr"
gds_file_name_2 <- ".gds"
### annotation file
dir_anno <- "/data/williamsjacr/UKB_WES_Full_Processed_Data/agds/"
anno_file_name_1 <- "Anno_chr"
anno_file_name_2 <- "_STAARpipeline.csv"

chr <- as.numeric(commandArgs(TRUE)[1])


###########################################################################
#           Main Function 
###########################################################################

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(readr)

### read annotation data
FunctionalAnnotation <- read_csv(paste0(dir_anno,"chr",chr,"/",anno_file_name_1,chr,anno_file_name_2),
                                 col_types=list(col_character(),col_double(),col_double(),col_double(),col_double(),
                                                col_double(),col_double(),col_double(),col_double(),col_double(),
                                                col_character(),col_character(),col_character(),col_double(),col_character(),
                                                col_character(),col_character(),col_character(),col_character(),col_double(),
                                                col_double(),col_character()))

dim(FunctionalAnnotation)

## rename colnames
colnames(FunctionalAnnotation)[2] <- "apc_conservation"
colnames(FunctionalAnnotation)[7] <- "apc_local_nucleotide_diversity"
colnames(FunctionalAnnotation)[9] <- "apc_protein_function"

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path, readonly = FALSE)

Anno.folder <- addfolder.gdsn(index.gdsn(genofile, "annotation/info"), "FunctionalAnnotation")
add.gdsn(Anno.folder, "FunctionalAnnotation", val=FunctionalAnnotation, compress="LZMA_ra", closezip=TRUE)

seqClose(genofile)