<!--
[RICE-ANNOTATION] Simulation_Study5
Purpose: Directory-level documentation for Simulation study pipelines (common + rare PRS on simulated phenotypes).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# Simulation_Study5 (simulation pipeline)

This folder runs one simulation scenario (“Simulation5” in the original output directory naming) end-to-end:
GWAS on simulated train data → common-variant PRSs (CT/LDpred2/lassosum2) → ensemble common PRS (RICE-CV)
→ rare-variant null models + STAAR gene-centric tests → rare-variant PRS (RICE-RV)
→ combined RICE model.

## How this maps to the paper
- **Manuscript**: Simulation Study methods and results (Figure 3).
- **Supplementary Figures**: Supp. Fig. 1–2.
- **Supplementary Data**: simulation sample sizes and related summary tables.

## Entry points
- `all.sh` : main driver (SLURM array) that orchestrates steps 1–7.
- `Common_Plus_Rare_PRS.R` : the final “combine” step called by `all.sh`.

## Subdirectories
- `GWAS_SummaryStatistics/` : GWAS sumstats from train data
- `CommonVariant_PRS/` : CT, LDpred2/lassosum2, and ensemble RICE-CV
- `RareVariant_Analysis/` : STAARpipeline preparation and gene-centric tests
- `RareVariant_PRS/` : rare-variant PRS construction

See the READMEs inside each subdirectory and the `[RICE-ANNOTATION]` headers for script-specific details.
