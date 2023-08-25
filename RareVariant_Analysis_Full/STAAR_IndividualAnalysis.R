rm(list=ls())
gc()

############## Input
## Null Model
obj_nullmodel <- get(load("/data/BB_Bioinformatics/ProjectData/UKB_WES_lipids/obj.STAAR.UKB.LDL.20211014.Rdata"))

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
output_path <- "/data/williamsjacr/UKB_WES_lipids/STAAR/Individual_Analysis/LDL/Results/"
# ## output file name
# output_file_name <- "UKBB_WES_200k_individual_analysis_LDL"

## output file name
output_file_name <- "UKBB_WES_200k_individual_analysis_LDL_Jake"


############## load source code
## load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)


#### Main Code
chr <- 1

# ### gds file
# gds.dir <- "/data/BB_Bioinformatics/ProjectData/UKB_WES_lipids/GDS/"
# gds.path <- paste0(gds.dir,"ukbb_wes_200k_chr",chr,".gds")
# genofile <- seqOpen(gds.path)

### gds file
gds.dir <- "/data/williamsjacr/UKB_WES_lipids/Data/gds/"
gds.path <- paste0(gds.dir,"all_chr",chr,".gds")
genofile <- seqOpen(gds.path)

position <- as.numeric(seqGetData(genofile, "position"))
start_loc <- min(position)
end_loc <- max(position)

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc < end_loc){
  results_individual_analysis <- Individual_Analysis(chr=chr, start_loc=start_loc, end_loc=end_loc, genofile=genofile, obj_nullmodel=obj_nullmodel, variant_type=variant_type, QC_label=QC_label, geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",chr,".Rdata"))

seqClose(genofile)



