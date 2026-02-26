<!--
[RICE-ANNOTATION] Simulation_Study8
Purpose: Directory-level documentation for Simulation study pipelines (common + rare PRS on simulated phenotypes).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# Simulation_Study8 (simulation pipeline)

End-to-end pipeline for the “Simulation8” scenario. This scenario is typically paired with phenotypes produced by
`Simulate_Data/Simulate_Data_NewProp_CausalRareVariants.R` (see `Simulate_Data/`), but uses a distinct simulated phenotype set under the
“Simulation8” output directory naming.

The folder structure matches `Simulation_Study5/`, `Simulation_Study6/`, and `Simulation_Study7/`.

## Entry points
- `all.sh` : main driver (SLURM array) that orchestrates the full pipeline.
- `Common_Plus_Rare_PRS.R` : final combine step.

See `[RICE-ANNOTATION]` headers for paper linkage and environment/path notes.
