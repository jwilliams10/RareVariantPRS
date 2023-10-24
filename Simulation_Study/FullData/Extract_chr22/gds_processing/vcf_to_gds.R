print(commandArgs(TRUE))
chr <- as.numeric(commandArgs(TRUE)[1])

library(SeqArray)

if(!file.exists(paste0("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds",chr,".gds"))){
    SeqArray::seqVCF2GDS("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_rare.vcf.bgz", out.fn = paste0("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds",chr,".gds"), header = NULL, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)
}