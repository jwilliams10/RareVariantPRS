<!--
[RICE-ANNOTATION] Simulation_Study6
Purpose: Directory-level documentation for Simulation study pipelines (common + rare PRS on simulated phenotypes).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# Simulation_Study6 (simulation pipeline)

This folder is structurally identical to `Simulation_Study5/` but runs the “Simulation6” scenario (distinct simulated phenotype set / architecture as produced in the original output directory naming).

## How this maps to the paper
- **Manuscript**: Simulation Study methods and results (Figure 3).
- **Supplementary Figures**: Supp. Fig. 1–2.
- **Supplementary Data**: simulation sample sizes and related summary tables.

## Entry points
- `all.sh` : main driver (SLURM array) that orchestrates steps 1–7.
- `Common_Plus_Rare_PRS.R` : the final “combine” step called by `all.sh`.

## Subdirectories
- `GWAS_SummaryStatistics/`
- `CommonVariant_PRS/`
- `RareVariant_Analysis/`
- `RareVariant_PRS/`

See the READMEs inside each subdirectory and the `[RICE-ANNOTATION]` headers for script-specific details.
