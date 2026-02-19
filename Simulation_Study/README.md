# Simulation_Study (UK Biobank WES chr22)

This directory contains the end-to-end pipeline used for the **chromosome 22 simulation study** in the  **RICE manuscript**. The code:
1) prepares chromosome 22 common and rare variant data from UKB WES,
2) simulates quantitative phenotypes under multiple causal/heritability scenarios,
3) creates train/tune/validation splits,
4) runs GWAS and common-variant PRS methods,
5) runs rare-variant association testing (STAARpipeline) and constructs rare-variant PRSs, and
6) combines common + rare PRSs into the final **RICE** predictor.

---

## How this maps to the paper

Primary references used throughout the inline annotations in this folder:

- **Manuscript:** 
  - Methods → *Simulation Study* (simulation design, causal proportions, negative selection/standardization)  
  - Results → *Simulation Study Results* (main simulation figure)

- **Supplementary data:**
  - Supplementary Data 1 / sheet **“S1 Sample Sizes; Sim. Study”** (train/tune/validation sample sizes)

- **Supplementary figures + note:**
  - **Supplementary Fig. 1–2** (additional simulation results + heritability checks)  
  - Supplementary Note → **“Ancestry Adjusted PRS”** (PC-based mean/variance adjustment implemented in downstream PRS scripts)

---

## Simulation design summary (paper defaults)

The simulation uses **unrelated** UKB WES participants and is based on **chromosome 22**. Causality is assigned to:
- randomly selected **common variants** (MAF > 0.01), and
- randomly selected **rare-variant sets** (gene-centric coding burdens on chr22).

Key scenario dimensions (see manuscript Methods and Supplementary Fig. 2):
- **Causal proportions** for both common variants and rare-variant sets: 0.01, 0.05, 0.20  
  - (Some scripts also include additional smaller levels, e.g., 0.001 and 0.0005.)
- **Negative selection / standardization (“Scaled: Yes/No”)** is implemented by whether genotype/burden matrices are standardized before applying effect sizes.
- **Within-set rare-variant causality**:
  - *All causal in set*: the burden for a causal set uses all rare variants in that set.
  - *Proportion causal in set*: a random proportion of variants within each causal set are treated as causal (see the `*_Prop_CausalRareVariants*` scripts).

Heritability settings used in the manuscript-aligned (“NewProp”) simulation scripts:
- **Common-variant component:** `h2_common = 0.05`
- **Rare-variant set component:** `h2_rare = 0.05 / 12 ≈ 0.00417`

> Note: `Simulate_Data/Simulate_Data.R` and the non-“NewProp” scripts use an older legacy setting (`h2_rare = 0.05/4`). The “NewProp” scripts match the manuscript text and Supplementary Fig. 2 description.

---

## Sample sizes (Supplementary Data 1)

From Supplementary Data → **S1 Sample Sizes; Sim. Study** (reported ancestries: AFR/AMR/EUR/SAS):

| Split | EUR | AFR | AMR | SAS |
|---|---:|---:|---:|---:|
| Train | 98,343 | – | – | – |
| Tune | 15,847 | 1,640 | 1,542 | 1,840 |
| Validation | 15,847 | 1,641 | 1,541 | 1,839 |

The split scripts in `Train_Tune_Validation_Split/` implement:
- **70% EUR training** (≈98k; used for the main simulation results), and
- a smaller **35% EUR training** variant (≈49k; used for sensitivity analyses).

The raw underlying dataset also contains a small EAS group; EAS is tracked in some scripts but is not included in the paper’s reported simulation performance tables/figures.

---

## Folder layout

- `FullData/`
  - `Extract_chr22/` – scripts to subset chr22 from UKB WES, create GDS/AGDS, add annotations/QC labels (data prep for simulation).
  - `G_star/` – functions and wrappers to compute rare-variant set burden matrices (G*) for gene-centric and sliding window set definitions.

- `Simulate_Data/` – phenotype simulation scripts (all-causal vs proportion-causal rare variants; legacy vs NewProp heritability ratio).

- `Train_Tune_Validation_Split/` – constructs train/tune/validation splits and writes per-split phenotype files used downstream.

- `GWAS_SummaryStatistics/` – runs PLINK2 GWAS in training data for each simulated phenotype replicate.

- `CommonVariant_PRS/` – common-variant PRS methods (C+T, LDpred2, lassosum2) and the common-variant ensemble (RICE-CV).

- `RareVariant_Analysis/` – fits STAAR null models and runs STAARpipeline association tests (gene-centric coding is enabled in `all.sh`).

- `RareVariant_PRS/` – constructs the rare-variant PRS (RICE-RV) using rare-variant set results.

- `Common_Plus_Rare_PRS.R` – final integration/evaluation script that produces **RICE** and computes performance (with ancestry-adjusted PRS standardization).

---

## Running the pipeline (high level)

This codebase is designed for an HPC/SLURM environment and includes many **hard-coded paths** under `/data/williamsjacr/...`.

Typical execution flow:
1) (One-time) Prepare chr22 data and G* burdens: `FullData/Extract_chr22/` and `FullData/G_star/`
2) Simulate phenotypes: `Simulate_Data/*.R` (run via the corresponding `*.sh` SLURM wrappers)
3) Create splits: `Train_Tune_Validation_Split/*.R`
4) Run the evaluation pipeline: `all.sh` (SLURM array; one task per simulated phenotype replicate index `i`)

See the per-file **[RICE-ANNOTATION]** headers for inputs/outputs and how each script links back to the manuscript and supplementary materials.
