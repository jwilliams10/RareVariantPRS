# =============================================================================
# [RICE-ANNOTATION] Simulate_Data/FullData/Extract_chr22/gds_processing/vcf_to_gds.R
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
print(commandArgs(TRUE))
chr <- as.numeric(commandArgs(TRUE)[1])

library(SeqArray)

if(!file.exists(paste0("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds",chr,".gds"))){
    SeqArray::seqVCF2GDS("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_rare.vcf.bgz", out.fn = paste0("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds",chr,".gds"), header = NULL, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)
}