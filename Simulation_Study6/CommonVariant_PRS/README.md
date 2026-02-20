<!--
[RICE-ANNOTATION] Simulation_Study6/CommonVariant_PRS
Purpose: Directory-level documentation for Simulation study pipelines (common + rare PRS on simulated phenotypes).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# CommonVariant_PRS

    Scripts to run common-variant PRS methods (CT, LDpred2, lassosum2) and the ensemble common-variant PRS (RICE-CV) within this simulation scenario.

    ## How this maps to the paper
    - Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.

    ## Files
    - `CT.R`: Runs clumping-and-thresholding (CT) PRS construction/evaluation for common variants.
- `LDPred_LASSOSum.R`: Runs LDpred2 and/or lassosum2 PRS construction/evaluation for common variants.
- `OneCommonPRS_All.R`: Builds the ensemble common-variant PRS (RICE-CV) by combining multiple PRS methods.
