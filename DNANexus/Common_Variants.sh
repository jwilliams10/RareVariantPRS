# =============================================================================
# [RICE-ANNOTATION] DNANexus/Common_Variants.sh
# Purpose: Shell wrapper to execute the paired R script on the target compute environment.
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
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; 
do
  dx run app-swiss-army-knife -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.fam -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bed -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bim -icmd="plink -bfile ukb24305_c${chrom}_b0_v1 --maf 0.01 --hwe 0.000001 --geno 0.02 --mind 0.05 --make-bed --out chr${chrom}_filtered_common" -y --destination UKB_PRS:JW/Clean_Data/ --instance-type mem1_ssd1_v2_x72
done

for chrom in 1 2 3 4 5; 
do
  dx run app-swiss-army-knife -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.fam -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bed -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bim -icmd="plink -bfile ukb24305_c${chrom}_b0_v1 --maf 0.01 --hwe 0.000001 --geno 0.02 --mind 0.05 --make-bed --out chr${chrom}_filtered_common" -y --destination UKB_PRS:JW/Clean_Data/ --instance-type mem1_hdd1_v2_x36
done