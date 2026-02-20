# =============================================================================
# [RICE-ANNOTATION] Imputed/reference_data.R
# Purpose: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
#
# Paper linkage:
#   - Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
#   - Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
#   - Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list = ls())

library(bigsnpr)

system(paste0("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/reference.txt --maf 0.01 --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference"))

if(file.exists("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.bk")){
  file.remove("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.bk")
}
if(file.exists("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.rds")){
  file.remove("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.rds")
}
snp_readBed("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference.bed",backingfile = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/all_chr_EUR_reference")