#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=10G

gunzip -c /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr22/ukbb_wes_200k_chr22.vcf.bgz > /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr22/ukbb_wes_200k_chr22.vcf