#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G

# =============================================================================
# [RICE-ANNOTATION] FullData/Extract_chr22/RareVariants.sh
# Purpose: Filters the chromosome-22 WES pVCF to the selected unrelated sample IDs to create a rare-variant VCF used for GDS/annotation and burden construction.
#
# Paper linkage:
#   - Manuscript; Methods → Simulation Study; Results → Simulation Study Results (Fig. 3).
#   - Supplementary Data: Supplementary Data 1 (sheet “S1 Sample Sizes; Sim. Study”).
#   - Supplementary Figures: Supp. Fig. 1–2; Supp. Note → “Ancestry Adjusted PRS”.
# =============================================================================

module load samtools

bcftools view -S /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_rare.txt --force-samples -o /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_rare.vcf.bgz /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr22/ukbb_wes_200k_chr22.vcf.bgz