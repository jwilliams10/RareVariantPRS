<!--
[RICE-ANNOTATION] Simulate_Data/FullData/Extract_chr22/gds_processing
Purpose: Directory-level documentation for Simulation data generation and chr22 WES data prep (used by Simulation_Study*).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# gds_processing (VCF → GDS/AGDS + annotations)

This folder converts the extracted chr22 rare-variant VCF into GDS/AGDS files and attaches analysis-ready fields (QC labels, functional annotations, and catalogs used by STAARpipeline).

## Typical workflow
- `all.sh` runs the full pipeline (often via SLURM job arrays):
  - `vcf_to_gds.R` : VCF → GDS
  - `Varinfo_gds.R` : add variant info tables
  - `Annotate.R` : attach functional annotation channels (used in weighting)
  - `Add_QC_label.R` : attach QC labels
  - `gds2agds.R` : shard into AGDS for parallel rare-variant analysis
  - `Association_Analysis_Prestep.R` : generate helper catalogs/metadata for downstream rare-variant steps
