<!--
[RICE-ANNOTATION] Imputed/CommonVariant_PRS
Purpose: Directory-level documentation for UKB imputed + WES analyses.

Paper linkage:
- Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
-->
# CommonVariant_PRS

Scripts to build common-variant PRSs and the ensemble common-variant model (**RICE-CV**) for the UKB imputed + WES pipeline.

## How this maps to the paper
- Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
- Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
- Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).

## Files
- `CT.R`: Runs clumping-and-thresholding (CT) PRS construction/evaluation for common variants.
- `CT.sh`: SLURM wrapper for a common-variant PRS method run (CT / LDpred2 / lassosum2 / ensemble).
- `CT_NewPCs.R`: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
- `CT_NewPCs.sh`: SLURM wrapper for a common-variant PRS method run (CT / LDpred2 / lassosum2 / ensemble).
- `Extract_Beta_All.R`: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
- `Extract_Beta_All.sh`: Shell wrapper to execute the paired R script on the target compute environment.
- `LDPred_LASSOSum.R`: Runs LDpred2 and/or lassosum2 PRS construction/evaluation for common variants.
- `LDPred_LASSOSum.sh`: SLURM wrapper for a common-variant PRS method run (CT / LDpred2 / lassosum2 / ensemble).
- `OneCommonPRS_All.R`: Builds the ensemble common-variant PRS (RICE-CV) by combining multiple PRS methods.
- `OneCommonPRS_All.sh`: SLURM wrapper for a common-variant PRS method run (CT / LDpred2 / lassosum2 / ensemble).
