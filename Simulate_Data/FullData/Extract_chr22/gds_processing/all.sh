#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=22
#SBATCH --mem-per-cpu=10G

# =============================================================================
# [RICE-ANNOTATION] FullData/Extract_chr22/gds_processing/all.sh
# Purpose: Driver shell script that runs the GDS processing pipeline (VCF→GDS→annotation/QC→AGDS) for chr22.
#
# Paper linkage:
#   - Manuscript; Methods → Simulation Study; Results → Simulation Study Results (Fig. 3).
#   - Supplementary Data: Supplementary Data 1 (sheet “S1 Sample Sizes; Sim. Study”).
#   - Supplementary Figures: Supp. Fig. 1–2; Supp. Note → “Ancestry Adjusted PRS”.
# =============================================================================

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/FullData/Extract_chr22/gds_processing/vcf_to_gds.R ${SLURM_ARRAY_TASK_ID} > vcf_to_gds"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/FullData/Extract_chr22/gds_processing/Add_QC_label.R > Add_QC_label.Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/FullData/Extract_chr22/gds_processing/Varinfo_gds.R ${SLURM_ARRAY_TASK_ID} > Varinfo_gds"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/FullData/Extract_chr22/gds_processing/Annotate.R ${SLURM_ARRAY_TASK_ID} > Annotate"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/FullData/Extract_chr22/gds_processing/gds2agds.R ${SLURM_ARRAY_TASK_ID} > gds2agds"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/FullData/Extract_chr22/gds_processing/Association_Analysis_Prestep.R > Association_Analysis_Prestep.Rout