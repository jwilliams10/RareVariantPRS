<!--
[RICE-ANNOTATION] README.md
Purpose: Directory-level documentation for Manuscript figures/tables and cross-dataset comparisons.

Paper linkage:
- Manuscript: Figure/table generation and cross-platform comparisons.
- Supplementary Figures/Data: see script-specific outputs and top-level README.
-->
# RareVariantPRS (RICE-CV / RICE-RV)

This repository contains the analysis pipelines and figure/table generation scripts used for the associated **Manuscript**. The code implements:
- **RICE-CV**: an ensemble common-variant PRS (combining CT, LDpred2, and lassosum2).
- **RICE-RV**: a rare-variant PRS built from gene-centric burden sets and penalized regression.
- **RICE**: the joint model combining RICE-CV and RICE-RV.

All scripts have been annotated with a standardized header tagged **`[RICE-ANNOTATION]`** that:
- states the file’s purpose,
- points readers to the relevant parts of the **Manuscript**, **Supplementary Data**, and **Supplementary Figures**, and
- notes where environment-specific paths (HPC/DNAnexus/AoU workbench) may need edits.

## Repository layout (high-level)

- `Imputed/`  
  UK Biobank **imputed + WES** analysis pipeline (common PRS, rare-variant tests/PRS, combined RICE).  
  Used for **Manuscript** Figures 4–5 and related **Supplementary Figures** / **Supplementary Data**.

- `DNANexus/`  
  UK Biobank **WGS** analysis pipeline intended to run on **DNAnexus**.  
  Used for UKB WGS results (context for **Manuscript** Figure 2) and WGS vs Imputed+WES comparisons (Figure 6), plus related **Supplementary Figures** / **Supplementary Data**.

- `AoU/`  
  **All of Us** (AoU) WGS analysis pipeline and portability evaluation.  
  Used for **Manuscript** Figures 7–8 and related **Supplementary Figures** / **Supplementary Data**.

- `Simulate_Data/`  
  Simulation data generation utilities (phenotype simulation + chr22 WES data preparation).  
  These scripts support the simulation pipelines below (Figure 3; Supplementary Figures 1–2).

- `Simulation_Study5/` – `Simulation_Study8/`  
  End-to-end simulation pipelines (GWAS → common PRS → rare-variant analysis/PRS → combined RICE) for different simulated architectures.  
  Each folder has an `all.sh` driver that orchestrates the run.

## Figure / table assembly scripts (top-level)

These scripts typically **read results produced by the pipelines above** and create the paper-ready plots/tables:
- `Fig2.R` (Manuscript Figure 2)
- `Simulation_Results.R` (Manuscript Figure 3; Supplementary Figure 1)
- `Sim_Characteristics.R` (Supplementary Figure 2)
- `WES_vs_WGS_Results*.R` (Manuscript Figure 6 and supporting Supplementary Figures)
- `Metrics_Table.R`, `LDSC_Results.R`, `Simulation_SampleSize.R` (Supplementary Data workbook population)
- `QQPlots.R` (association diagnostic QQ plots used in Supplementary Figures)

## Notes on execution

Many scripts contain absolute paths for the original compute environments (HPC scratch, DNAnexus project mounts, AoU Researcher Workbench).
The annotations and per-folder READMEs call this out and describe expected inputs/outputs so paths can be updated for a new setup.
