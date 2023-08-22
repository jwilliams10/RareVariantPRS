#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=10G

module load vcftools

vcftools --gzvcf /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_train.vcf.bgz --freq2 --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_train --max-alleles 2
vcftools --gzvcf /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_tune.vcf.bgz --freq2 --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_tune --max-alleles 2
vcftools --gzvcf /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_validation.vcf.bgz --freq2 --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_validation --max-alleles 2