<!--
[RICE-ANNOTATION] Simulation_Study6/RareVariant_Analysis
Purpose: Directory-level documentation for Simulation study pipelines (common + rare PRS on simulated phenotypes).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# RareVariant_Analysis

Scripts to fit rare-variant null models and run STAAR gene-centric tests within this simulation scenario.

## How this maps to the paper
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.

## Files
- `NullModel.R`: Fits null models required for gene-centric rare variant tests (STAARpipeline / burden tests).
- `STAAR_GeneCentric_Coding.R`: Runs gene-centric rare-variant association testing using STAAR/STAARpipeline outputs.
