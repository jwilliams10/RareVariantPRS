#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G

module load samtools

bcftools view -S /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/sampleids_rare.txt --force-samples -o /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_rare.vcf.bgz /data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr22/ukbb_wes_200k_chr22.vcf.bgz