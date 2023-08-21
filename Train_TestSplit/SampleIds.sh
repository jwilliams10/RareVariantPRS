#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=10G

module load samtools

bcftools query -l /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr21/ukbb_wes_200k_chr21_test.vcf.bgz > /data/williamsjacr/UKB_WES_lipids/Data/attempt_21_tune.txt