# =============================================================================
# [RICE-ANNOTATION] Imputed/DataSpecific_PCs.R
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

system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --maf 0.05 --indep-pairwise 50 5 0.8 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/train_pc_dat")
system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/train.txt --extract /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/train_pc_dat.prune.in --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/train_pc_dat")
system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/train_pc_dat --pca 10 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/train_pc_dat --threads 20")


system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --maf 0.05 --indep-pairwise 50 5 0.8 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/tune_pc_dat")
system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/tune.txt --extract /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/tune_pc_dat.prune.in --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/tune_pc_dat")
system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/tune_pc_dat --pca 10 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/tune_pc_dat --threads 20")

system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --maf 0.05 --indep-pairwise 50 5 0.8 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/validation_pc_dat")
system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --keep /data/williamsjacr/UKB_WES_Phenotypes/validation.txt --extract /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/validation_pc_dat.prune.in --make-bed --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/validation_pc_dat")
system("/data/williamsjacr/software/plink2 --bfile /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/validation_pc_dat --pca 10 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/validation_pc_dat --threads 20")