# =============================================================================
# [RICE-ANNOTATION] AoU/JointPRS_PROSPER_MUSSEL_SNPINFO.R
# Purpose: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
#
# Paper linkage:
#   - Manuscript: Results -> All of Us Results (Fig. 7) and AoU-trained PRS evaluated on UKB (Fig. 8).
#   - Supplementary Figures: Supp. Fig. 17–21 (AoU association diagnostics, PRS performance, and portability).
#   - Supplementary Data: key AoU cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list = ls())
library(data.table)

snpinfo_mult_1kg_hm3 <- read.delim("/data/williamsjacr/PRSCSx_LD/snpinfo_mult_1kg_hm3")
ref_bim_PROSPER <- fread("/data/williamsjacr/MUSSEL_PROSPER_LDs/PROSPER/ref_bim.txt")
ref_bim_PROSPER <- as.data.frame(ref_bim_PROSPER)
ref_bim_MUSSEL <- fread("/data/williamsjacr/MUSSEL_PROSPER_LDs/ref_bim.txt")
ref_bim_MUSSEL <- as.data.frame(ref_bim_MUSSEL)

SNP_Info_JointPRS_PROSPER_MUSSEL <- data.frame(SNP = c(ref_bim_MUSSEL$V2,ref_bim_PROSPER$V2,snpinfo_mult_1kg_hm3$SNP),
              CHR = c(ref_bim_MUSSEL$V1,ref_bim_PROSPER$V1,snpinfo_mult_1kg_hm3$CHR),
              BP = c(ref_bim_MUSSEL$V4,ref_bim_PROSPER$V4,snpinfo_mult_1kg_hm3$BP),
              A1 = c(ref_bim_MUSSEL$V5,ref_bim_PROSPER$V5,snpinfo_mult_1kg_hm3$A1),
              A2 = c(ref_bim_MUSSEL$V6,ref_bim_PROSPER$V6,snpinfo_mult_1kg_hm3$A2))

SNP_Info_JointPRS_PROSPER_MUSSEL <- unique(SNP_Info_JointPRS_PROSPER_MUSSEL)
SNP_Info_JointPRS_PROSPER_MUSSEL <- SNP_Info_JointPRS_PROSPER_MUSSEL[!duplicated(SNP_Info_JointPRS_PROSPER_MUSSEL$SNP),]

write.csv(SNP_Info_JointPRS_PROSPER_MUSSEL,file = "/data/williamsjacr/MUSSEL_PROSPER_LDs/SNP_Info_JointPRS_PROSPER_MUSSEL.csv",row.names = FALSE)