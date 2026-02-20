<!--
[RICE-ANNOTATION] Simulate_Data/FullData
Purpose: Directory-level documentation for Simulation data generation and chr22 WES data prep (used by Simulation_Study*).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# FullData (simulation inputs derived from UKB WES chr22)

This directory holds one-time preprocessing steps to extract and format UKB WES chromosome 22 data for the simulation study.

## Subdirectories
- `Extract_chr22/`  
  Select unrelated samples, build chr22 common-variant PLINK files, subset rare variants, and convert rare-variant VCF → GDS/AGDS (for STAARpipeline and burden construction).
- `G_star/`  
  Constructs gene-centric coding burden-score matrices (“G*”) used by the rare-variant PRS components in the simulation pipelines.

See per-file `[RICE-ANNOTATION]` headers and the READMEs inside each subdirectory for the recommended run order.
