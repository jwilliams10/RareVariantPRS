<!--
[RICE-ANNOTATION] AoU
Purpose: Directory-level documentation for All of Us (AoU) WGS analyses.

Paper linkage:
- Manuscript: Results -> All of Us Results (Fig. 7) and AoU-trained PRS evaluated on UKB (Fig. 8).
- Supplementary Figures: Supp. Fig. 17–21 (AoU association diagnostics, PRS performance, and portability).
- Supplementary Data: key AoU cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
-->
# AoU (All of Us) analyses

This folder contains scripts and notebooks used for the **All of Us (AoU) WGS** analyses described in the **Manuscript** (Results → *All of Us Results* and *Portability*).

## How this maps to the paper
- **Manuscript**: AoU results (Figure 7) and AoU-trained PRS evaluated in UKB (Figure 8).
- **Supplementary Figures**: Supp. Fig. 17–21 (AoU association diagnostics, PRS performance, and portability).
- **Supplementary Data**: cohort summaries and diagnostics (e.g., sample sizes, genomic control metrics, variant counts, runtime tables).

## Key entry points
- `Submission_Script.ipynb`  
  AoU Researcher Workbench notebook for job setup/submission (project- and workspace-specific settings live here).
- `Figures.ipynb`  
  Notebook used to compile AoU outputs into paper-ready plots (see inline `[RICE-ANNOTATION]` cell).
- `Overall_Results.R`  
  Aggregates AoU PRS performance results across traits/ancestries.
- `LDSC.R` / `LDSC.sh`  
  Runs LD score regression (or wraps LDSC calls) for AoU summary statistics.
- `JointPRS_Sumstats.R` / `JointPRS_Sumstats.sh` + `JointPRS_*.sh`  
  Prepares/feeds summary statistics to JointPRS; trait-specific wrappers are provided.
- `RICECV_AoU_CrossPlatform.*`, `RICERV_AoU_CrossPlatform.*`, `RICE_AoU_CrossPlatform.*`  
  Cross-platform evaluation utilities (AoU ↔ UKB comparisons).
- `AoU_CV_UKB_RV.*`  
  Runs analyses comparing common-variant AoU PRS with UKB rare-variant components.

## Conventions
- Many `.R` scripts have paired `.sh` wrappers for running on a scheduler.
- Paths are environment-specific (AoU workspace buckets / local mounts); update as needed.

See the per-file `[RICE-ANNOTATION]` headers for the intended role of each script.
