#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=0-41
#SBATCH --mem-per-cpu=10G

module load samtools

INPUT_PATH=/data/BB_Bioinformatics/ProjectData/UKB/exome_seq/OQFE_VCF
OUTPUT_PATH=/data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr9
CHR=9

## step1: multiallelic caller
bcftools norm -Oz -m - $INPUT_PATH/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_1.vcf.gz
## step2: calculate HWE
bcftools +fill-tags ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_1.vcf.gz -Oz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_2.vcf.gz -- -t HWE
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_1.vcf.gz
## step3: filter by HWE<1E-15
bcftools filter -Oz -e 'HWE<1E-15' --set-GTs . ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_2.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_3.vcf.gz
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_2.vcf.gz
## step4: filter by DP<7
bcftools filter -Oz -e 'DP<7' --set-GTs . ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_3.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_4.vcf.gz
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_3.vcf.gz
## step5: genotype filtering
bcftools filter -Oz -e '(GT="het" & AD[:0]/(AD[:0]+AD[:1])<0.15) | (GT="het" & AD[:1]/(AD[:0]+AD[:1])<0.15) | (GT="het" & GQ<20) | (GT="het" & binom(AD)<1E-3)' --set-GTs . ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_4.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_5.vcf.gz
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_4.vcf.gz
## step6: calculate missingness
bcftools +fill-tags -Oz ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_5.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_6.vcf.gz -- -t F_MISSING
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_5.vcf.gz
## step7: remove variants with missing>50%
bcftools filter -e 'F_MISSING>0.5' -Oz ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_6.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_7.vcf.gz 
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_6.vcf.gz
## step8: remove FORMAT
bcftools annotate -x FORMAT/DP,FORMAT/AD,FORMAT/GQ,FORMAT/PL,FORMAT/RNC -Oz ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_7.vcf.gz -o ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_8.vcf.gz 
rm ${OUTPUT_PATH}/ukb23156_c${CHR}_b${SLURM_ARRAY_TASK_ID}_v1_7.vcf.gz