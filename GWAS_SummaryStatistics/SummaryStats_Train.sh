#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=100G

/data/williamsjacr/software/plink2 --pfile /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}_common_train \
                                   --pheno /data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt \
                                   --pheno-name LDLadj.norm \
                                   --linear \
                                   --covar /data/williamsjacr/UKB_WES_lipids/Data/phenotypes/LDL_Train.txt \
                                   --covar-name age, age2, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
                                   --vif 999 \
                                   --out /data/williamsjacr/UKB_WES_lipids/Data/GWAS_Summary_Statistics/LDL/LDL_adj_norm__chr${SLURM_ARRAY_TASK_ID}_train