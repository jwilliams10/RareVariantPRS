<!--
[RICE-ANNOTATION] Simulate_Data/FullData/G_star
Purpose: Directory-level documentation for Simulation data generation and chr22 WES data prep (used by Simulation_Study*).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# G_star (burden score construction for simulation)

This folder constructs gene-centric **coding** burden-score matrices (“G” / “G*”) from the chr22 rare-variant data prepared in `FullData/Extract_chr22/`.

## Key scripts
- `GSTAR_Extract_GeneCentric_Coding.R` / `.sh`  
  Main extractor for gene-centric coding sets (produces burden score inputs used in downstream simulation pipelines).
- `RareVariants_Extract_GeneCentric_Coding.R`  
  Helper extractor/formatter for coding rare variants prior to burden-score construction.
- `g_star_gene_centric_coding.R`  
  Utilities for generating/validating the gene-centric coding burden matrices.