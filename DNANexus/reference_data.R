# =============================================================================
# [RICE-ANNOTATION] DNANexus/reference_data.R
# Purpose: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
#
# Paper linkage:
#   - Manuscript: Results -> UKB WGS Results (Fig. 2 context; WGS analyses) and WGS vs Imputed+WES comparison (Fig. 6).
#   - Supplementary Figures: Supp. Fig. 5–6 (WGS association diagnostics) and Supp. Fig. 12–16 (PRS performance + coding/noncoding comparisons).
#   - Supplementary Data: key UKB WGS cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list = ls())
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/reference_data.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/reference_data.sh  -icmd="bash reference_data.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Data/ --instance-type mem1_ssd1_v2_x36

if(!("remotes" %in% rownames(installed.packages()))){
  install.packages("remotes",quiet = TRUE)
}

if(!("bigsnpr" %in% rownames(installed.packages()))){
  install.packages("bigsnpr",quiet = TRUE)
}

library(bigsnpr)

for(i in 1:22){
  system(paste0("plink2 --bfile Clean_Data/chr",i,"_filtered_common --keep reference.txt --make-bed --out chr",i,"_filtered_common_reference"))
}

system(paste0("plink2 --bfile Clean_Data/all_chr --keep reference.txt --make-bed --out all_chr_reference"))

bigsnpr::snp_readBed("all_chr_reference.bed",backingfile = "all_chr_reference")

system("rm -r Clean_Data/")
file.remove("reference.txt")