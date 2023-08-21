#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=10G

module load samtools

bcftools query -l /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1.vcf.bgz > /data/williamsjacr/UKB_WES_lipids/Data/sampleids.txt