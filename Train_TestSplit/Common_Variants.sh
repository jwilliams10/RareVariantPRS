#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=100G

/data/williamsjacr/software/plink2 --vcf /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1.vcf.bgz --vcf-half-call m dosage=DS --make-pgen --maf 0.01 --out /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr1/ukbb_wes_200k_chr1_common