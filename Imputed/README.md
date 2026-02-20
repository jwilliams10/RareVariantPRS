<!--
[RICE-ANNOTATION] Imputed
Purpose: Directory-level documentation for UKB imputed + WES analyses.

Paper linkage:
- Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
-->
# Imputed (UK Biobank imputed + WES pipeline)

This folder contains the UK Biobank **imputed genotype + WES** pipeline used to construct common-variant PRSs, rare-variant PRSs, and the combined RICE model for the analyses reported in the **Manuscript**.

## How this maps to the paper
- **Manuscript**: UKB Imputed + WES results (Figures 4–5).
- **Supplementary Figures**: association diagnostics and PRS performance / sensitivity (e.g., Supp. Fig. 3–4, 7–11).
- **Supplementary Data**: cohort summaries and diagnostics (e.g., sample sizes, genomic control metrics, variant counts, runtime tables).

## Subdirectories
- `CommonVariant_PRS/`  
  CT / LDpred2 / lassosum2 runs and the ensemble **RICE-CV** construction.
- `RareVariant_Analysis/`  
  Null models and STAARpipeline gene-centric rare variant association testing (coding masks, long masks).
- `RareVariant_PRS/`  
  Construction of **RICE-RV** and rare-variant sensitivity analyses.

## Top-level scripts in this folder
- `Common_Plus_Rare_PRS.*` : combine common + rare components into joint RICE
- `DataSpecific_PCs.*` : compute dataset-specific PCs used by “*_NewPCs” PRS runs
- `LDSC.*` : LD score regression / genomic inflation diagnostics
- `QQPlots_CV.R` : QQ plots for common-variant diagnostics
- `Overall_Results_*` : aggregate continuous/binary trait performance results
- `SampleSizes.R`, `Time_Table.R` : produce Supplementary Data tables
- `NRI_Tables*.R` : net reclassification improvement summaries (used in Results / Supplementary)

See per-file `[RICE-ANNOTATION]` headers for details (inputs/outputs and paper linkage).
