<!--
[RICE-ANNOTATION] Simulate_Data/FullData/Extract_chr22
Purpose: Directory-level documentation for Simulation data generation and chr22 WES data prep (used by Simulation_Study*).

Paper linkage:
- Manuscript: Methods -> Simulation Study; Results -> Simulation Study Results (Fig. 3).
- Supplementary Figures: Supp. Fig. 1–2 (simulation performance and design characteristics).
- Supplementary Data: simulation sample sizes and related summaries.
-->
# Extract_chr22 (UK Biobank WES chr22 data prep)

This folder contains the scripts used to prepare **chromosome 22** data from **UKB WES** for the simulation study described in the attached manuscript (Methods → *Simulation Study*).

The overall goal is to:
1) define an unrelated analysis sample,
2) subset chr22 common variants (for GWAS/PRS),
3) subset chr22 rare variants (for STAARpipeline / burden construction),
4) convert rare-variant VCF → GDS/AGDS and attach annotations/QC labels.

---

## Typical one-time workflow

### 1) Select unrelated individuals and build chr22 common-variant PLINK dataset
Run:
- `CommonVariants_SampleIDs.R`

Outputs (paths are hard-coded in the script):
- chr22 SNP list (HM3 subset)
- sample ID files for common vs rare variant extraction
- `/.../chr22_filtered_common.{bed,bim,fam}` (PLINK files used later for GWAS + common PRS)

### 2) Extract chr22 rare variants for the selected sample
Run (SLURM):
- `RareVariants.sh`

Output:
- `/.../chr22_filtered_rare.vcf.bgz`

### 3) Process the rare-variant VCF into analysis-ready GDS/AGDS
Run:
- `gds_processing/all.sh`

This wrapper drives:
- `gds_processing/vcf_to_gds.R` (VCF → GDS)
- `gds_processing/Varinfo_gds.R` + `gds_processing/Association_Analysis_Prestep.R` (build job catalogs / metadata used by STAARpipeline)
- `gds_processing/Annotate.R` (add functional annotation channels; used for STAAR weights)
- `gds_processing/Add_QC_label.R` (variant QC label field)
- `gds_processing/gds2agds.R` (split into AGDS shards for parallel analysis)

## Link to paper artifacts

Inline annotations in these scripts point to:
- Manuscript (Methods → Simulation Study)
- Supplementary Figures (Supplementary Fig. 2 describes “Scaled: Yes/No” and rare-variant set causality scenarios)
- Supplementary Data (Supplementary Data 1 sample sizes used later in the simulation splits)
