#!/bin/bash --login
# =============================================================================
# [RICE-ANNOTATION] AoU/JointPRS_LDL.sh
# Purpose: SLURM wrapper for running JointPRS for a specific trait or step.
#
# Paper linkage:
#   - Manuscript: Results -> All of Us Results (Fig. 7) and AoU-trained PRS evaluated on UKB (Fig. 8).
#   - Supplementary Figures: Supp. Fig. 17–21 (AoU association diagnostics, PRS performance, and portability).
#   - Supplementary Data: key AoU cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=30G

# module purge
module load python/3.8
source /usr/local/Anaconda/envs/py3.8/etc/profile.d/conda.sh && source /usr/local/Anaconda/envs/py3.8/etc/profile.d/mamba.sh

mamba activate JointPRS

reference_path=/data/williamsjacr/PRSCSx_LD; type=1KG
outcome_path=/data/williamsjacr/AoU_JointPRS
chr=${SLURM_ARRAY_TASK_ID}

pop1=EUR; pop2=AMR; pop3=AFR
r1=1; r2=1; r3=1
sst1=/data/williamsjacr/AoU_JointPRS/GWASSumStats/EUR_LDL_GWAS_SumStats_Cleaned.txt; sst2=/data/williamsjacr/AoU_JointPRS/GWASSumStats/AMR_LDL_GWAS_SumStats_Cleaned.txt; sst3=/data/williamsjacr/AoU_JointPRS/GWASSumStats/AFR_LDL_GWAS_SumStats_Cleaned.txt
sample_size1=42844; sample_size2=7824; sample_size3=11946

python /data/williamsjacr/software/JointPRS/JointPRS.py \
--ref_dir=${reference_path} \
--bim_prefix=/data/williamsjacr/AoU_JointPRS/GWASSumStats/EUR_LDL \
--pop=${pop1},${pop2},${pop3} \
--rho_cons=${r1},${r2},${r3} \
--sst_file=${sst1},${sst2},${sst3} \
--n_gwas=${sample_size1},${sample_size2},${sample_size3} \
--chrom=${chr} \
--out_dir=${outcome_path} \
--meta=TRUE \
--out_name=JointPRS_LDL