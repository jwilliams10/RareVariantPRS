#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G

# module purge
module load regenie/3.0.3

regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --phenoColList HDL,LDL,TC,logTG,BMI,Height --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --qt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_continuous &
wait
regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --phenoColList HDL,LDL,TC,logTG,BMI,Height --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_continuous_pred.list --qt --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous