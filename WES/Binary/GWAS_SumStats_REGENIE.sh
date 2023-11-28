#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G

# module purge
module load regenie/3.0.3

regenie --step 1 --bed /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --phenoColList Asthma,CAD,T2D,Breast,Prostate --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --bt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step1
regenie --step 2 --bed /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --phenoColList Asthma,CAD,T2D,Breast,Prostate --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/regenie_step1_pred.list --bt --firth --approx --pThresh 0.05 --bsize 400 --split --out /data/williamsjacr/UKB_WES_Phenotypes/Binary/GWAS_Summary_Statistics/sumstats