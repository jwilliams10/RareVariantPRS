#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=100G

/data/williamsjacr/software/plink2 --vcf /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}.vcf.bgz dosage=DS --vcf-half-call m --make-pgen --maf 0.01 --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common
/data/williamsjacr/software/plink2 --pfile /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common --keep /data/williamsjacr/UKB_WES_lipids/Data/train.txt --make-pgen --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common_train
/data/williamsjacr/software/plink2 --pfile /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common --keep /data/williamsjacr/UKB_WES_lipids/Data/tune.txt --make-pgen --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common_tune
/data/williamsjacr/software/plink2 --pfile /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common --keep /data/williamsjacr/UKB_WES_lipids/Data/validation.txt --make-pgen --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common_validation
