<!--
[RICE-ANNOTATION] Simulation_Study7
Purpose: Directory-level documentation for Simulation study pipelines (common + rare PRS on simulated phenotypes).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# Simulation_Study7 (simulation pipeline)

End-to-end pipeline for the “Simulation7” scenario. This scenario is typically paired with phenotypes produced by
`Simulate_Data_NewProp_CausalRareVariants.R` (see `Simulate_Data/`), i.e., architectures where only a proportion of rare variants
within a causal set are causal.

The folder structure matches `Simulation_Study5/6/8`.

## Entry points
- `all.sh` : main driver (SLURM array) that orchestrates the full pipeline.
- `Common_Plus_Rare_PRS.R` : final combine step.

See `[RICE-ANNOTATION]` headers for paper linkage and environment/path notes.
