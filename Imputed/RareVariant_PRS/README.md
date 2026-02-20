<!--
[RICE-ANNOTATION] Imputed/RareVariant_PRS
Purpose: Directory-level documentation for UKB imputed + WES analyses.

Paper linkage:
- Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
-->
# RareVariant_PRS

    Scripts to construct the rare-variant PRS component (**RICE-RV**) and related sensitivity analyses for the UKB imputed + WES pipeline.

    ## How this maps to the paper
    - Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).

    ## Files
    - `Lipids_BestGenes_Results.R`: Aggregates results across traits/ancestries and generates summary outputs.
- `RICE_RV_Lipids_Best.R`: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
- `RICE_RV_Lipids_Best.sh`: Shell wrapper to execute the paired R script on the target compute environment.
- `RICE_RV_Lipids_Gene.R`: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
- `RICE_RV_Lipids_Gene.sh`: Shell wrapper to execute the paired R script on the target compute environment.
- `Sensitivity_Analysis_RICE_RV.R`: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
- `Sensitivity_Analysis_RICE_RV.sh`: Shell wrapper to execute the paired R script on the target compute environment.
- `Sensitivity_Analysis_Results.R`: Aggregates results across traits/ancestries and generates summary outputs.
- `Single_RareVariant_PRS_All.R`: Constructs the rare-variant PRS component (RICE-RV) from burden scores and selected rare-variant sets.
- `Single_RareVariant_PRS_All.sh`: Shell wrapper to execute the paired R script on the target compute environment.
