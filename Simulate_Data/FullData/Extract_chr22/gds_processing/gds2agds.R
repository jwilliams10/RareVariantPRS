# =============================================================================
# [RICE-ANNOTATION] Simulate_Data/FullData/Extract_chr22/gds_processing/gds2agds.R
# Purpose: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
#
# Paper linkage:
#   - Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
#   - Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
#   - Supplementary Data: simulation sample sizes and related summaries.
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### gds file
dir_geno <- "/data/williamsjacr/UKB_WES_lipids/Data/gds/"
### annotation file
dir_anno <- "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/"
anno_file_name_1 <- "Anno_chr_train"
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
gds.path <- "/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds"
genofile <- seqOpen(gds.path, readonly = FALSE)

Anno.folder <- addfolder.gdsn(index.gdsn(genofile, "annotation/info"), "FunctionalAnnotation")
add.gdsn(Anno.folder, "FunctionalAnnotation", val=FunctionalAnnotation, compress="LZMA_ra", closezip=TRUE)

seqClose(genofile)

