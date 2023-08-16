#!/bin/bash
#SBATCH -J UKBBC
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-144:00
#SBATCH --array=5-7
#SBATCH -p test
#SBATCH --mem=10000
#SBATCH -o /n/holyscratch01/xlin/xihao_zilin/UKBB/pVCF/hostname_%j.out  
#SBATCH -e /n/holyscratch01/xlin/xihao_zilin/UKBB/pVCF/hostname_%j.err 
#SBATCH --mail-type=NONE

OUTPUT_PATH=/n/holyscratch01/xlin/xihao_zilin/UKBB/pVCF

/n/holystore01/LABS/xlin/Lab/xihao_zilin/UKB_WES_lipids/QC/bcftools/bin/bcftools concat ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/ukb23156_c${SLURM_ARRAY_TASK_ID}_b*_v1_8.vcf.gz -Oz -o ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}.vcf.bgz
