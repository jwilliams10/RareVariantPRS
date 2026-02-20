<!--
[RICE-ANNOTATION] DNANexus
Purpose: Directory-level documentation for UK Biobank (UKB) WGS analyses executed on DNAnexus.

Paper linkage:
- Manuscript: Results -> UKB WGS Results (Fig. 2 context; WGS analyses) and WGS vs Imputed+WES comparison (Fig. 6).
- Supplementary Figures: Supp. Fig. 5–6 (WGS association diagnostics) and Supp. Fig. 12–16 (PRS performance + coding/noncoding comparisons).
- Supplementary Data: key UKB WGS cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
-->
# DNANexus (UK Biobank WGS pipeline)

This folder contains the UK Biobank **WGS** analysis pipeline intended to run on **DNAnexus**. The scripts implement the end-to-end workflow:
GWAS summary statistics → common-variant PRS methods → RICE-CV (ensemble) → rare-variant analysis/PRS (RICE-RV) → combined RICE → result extraction/plots.

## How this maps to the paper
- **Manuscript**: UKB WGS analyses (context for Figure 2) and WGS vs Imputed+WES comparisons (Figure 6).
- **Supplementary Figures**: WGS diagnostics/performance (e.g., Supp. Fig. 5–6, 12–16).
- **Supplementary Data**: cohort summaries and diagnostics (e.g., sample sizes, genomic control metrics, variant counts, runtime tables).

## Pipeline structure (main components)
**GWAS / diagnostics**
- `GWAS_SumStats_Continuous.R`, `GWAS_SumStats_Binary.R` (generate WGS GWAS summary statistics)
- `LDSC.R` / `LDSC.sh` (LD score regression / genomic inflation diagnostics)
- `QQPlots_*` (QQ plots for association diagnostics; used in Supplementary Figures)

**Common-variant PRS (RICE-CV inputs)**
- `CT*.R` / `CT*.sh` (clumping & thresholding)
- `LDPred_LASSOSum*.R` / `.sh` (LDpred2 + lassosum2)
- `OneCommonPRS_All*.R` / `.sh` (ensemble common PRS construction)

**Rare-variant analysis and PRS (RICE-RV)**
- `NullModel*.R` / `.sh` (null models for gene-centric rare variant testing)
- `G_Extraction_*.R` / `.sh` (burden score extraction for coding/noncoding sets)
- `RV_Analysis_Summary*.R` / `.sh` (collate rare-variant association results)

**Joint model + reporting**
- `Common_Plus_Rare_PRS*.R` / `.sh` (combine RICE-CV + RICE-RV into joint RICE model)
- `Extract_Betas_All*.R` / `.sh`, `Extract_Results.R`, `Overall_Results*.R` (compile final metrics)
- `Coding_vs_Noncoding.R` (coding vs noncoding comparisons used in Supplementary Figures)

## Notes
- Many scripts are duplicated with `_Binary` suffix for binary traits.
- DNAnexus-specific mounts/paths are present in the scripts; update to match your DNAnexus project setup.

See the per-file `[RICE-ANNOTATION]` headers for script-specific intent and paper linkage.
