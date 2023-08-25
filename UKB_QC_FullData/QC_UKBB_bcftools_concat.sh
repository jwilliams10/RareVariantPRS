#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=19
#SBATCH --mem-per-cpu=10G

module load samtools

OUTPUT_PATH=/data/williamsjacr/UKB_WES_lipids/Data/pVCF

ls ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/ukb23156_c${SLURM_ARRAY_TASK_ID}_b*_v1_8.vcf.gz | sort -V > ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/sort.temp
bcftools concat --file-list ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/sort.temp -Oz -o ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/ukbb_wes_200k_chr${SLURM_ARRAY_TASK_ID}.vcf.bgz
rm ${OUTPUT_PATH}/chr${SLURM_ARRAY_TASK_ID}/sort.temp
