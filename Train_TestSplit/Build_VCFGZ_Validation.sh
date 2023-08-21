#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=10G

module load samtools

bcftools view -S /data/williamsjacr/UKB_WES_lipids/Data/validation.txt --force-samples -o /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_validation.vcf.bgz /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}.vcf.bgz