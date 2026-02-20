<!--
[RICE-ANNOTATION] Simulate_Data
Purpose: Directory-level documentation for Simulation data generation and chr22 WES data prep (used by Simulation_Study*).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# Simulate_Data (simulation data generation and chr22 prep)

This folder contains utilities used to:
1) generate simulated phenotypes under the genetic architectures described in the **Manuscript** (Methods → *Simulation Study*), and
2) prepare chr22 UKB WES data inputs used by those simulations.

## How this maps to the paper
- **Manuscript**: Simulation Study methods and results (Figure 3).
- **Supplementary Figures**: Supp. Fig. 1–2 (simulation performance + design characteristics).
- **Supplementary Data**: simulation sample sizes and related summary tables.

## Key scripts
- `Simulate_Data_NewProp.R` / `.sh`  
  Simulates phenotypes where all rare variants within a “causal set” contribute (architecture variant 1).
- `Simulate_Data_NewProp_CausalRareVariants.R` / `.sh`  
  Simulates phenotypes where only a *proportion* of rare variants within a causal set contribute (architecture variant 2).
- `Train_Tune_Validate_Split_NewProp*.R`  
  Creates train/tune/validation splits and per-split phenotype files for downstream simulation pipelines.

## Data preparation
- `FullData/Extract_chr22/` prepares chr22 common + rare variant inputs and converts rare-variant VCF → GDS/AGDS.
- `FullData/G_star/` constructs burden-score (“G matrix”) inputs for gene-centric coding sets.

