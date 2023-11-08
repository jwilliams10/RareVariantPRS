#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=100G

/data/williamsjacr/software/plink2 --vcf /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}.vcf.bgz dosage=DS --vcf-half-call m --make-bed --maf 0.01 --out /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common