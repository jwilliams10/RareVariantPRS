<!--
[RICE-ANNOTATION] Imputed/RareVariant_Analysis
Purpose: Directory-level documentation for UKB imputed + WES analyses.

Paper linkage:
- Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
-->
# RareVariant_Analysis

Scripts to fit null models and run gene-centric rare-variant association testing (STAAR/STAARpipeline) for the UKB imputed + WES pipeline.

## How this maps to the paper
- Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).

## Files
- `Modified_Summary_Script.R`: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
- `NullModels.R`: Fits null models required for gene-centric rare variant tests (STAARpipeline / burden tests).
- `STAAR_GeneCentric_Coding.R`: Runs gene-centric rare-variant association testing using STAAR/STAARpipeline outputs.
- `STAAR_GeneCentric_Coding.sh`: SLURM wrapper to run gene-centric rare-variant association tests (STAARpipeline / summaries).
- `STAAR_GeneCentric_Coding_LongMasks.R`: Runs gene-centric rare-variant association testing using STAAR/STAARpipeline outputs.
- `STAAR_GeneCentric_Coding_LongMasks.sh`: SLURM wrapper to run gene-centric rare-variant association tests (STAARpipeline / summaries).
