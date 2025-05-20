#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --array=1-5
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G

# module purge
module load regenie/3.0.3

if [ ${SLURM_ARRAY_TASK_ID} = 1 ]
then
       trait=Asthma
       regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --bt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait} &
       wait
       regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait}_pred.list --bt --firth --approx --pThresh 0.05 --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_binary_${trait}
elif [ ${SLURM_ARRAY_TASK_ID} = 2 ]
then
       trait=CAD
       regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --bt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait} &
       wait
       regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait}_pred.list --bt --firth --approx --pThresh 0.05 --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_binary_${trait}
elif [ ${SLURM_ARRAY_TASK_ID} = 3 ]
then
       trait=T2D
       regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --bt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait} &
       wait
       regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait}_pred.list --bt --firth --approx --pThresh 0.05 --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_binary_${trait}
elif [ ${SLURM_ARRAY_TASK_ID} = 4 ]
then 
       trait=Breast
       regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --bt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait} &
       wait
       regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait}_pred.list --bt --firth --approx --pThresh 0.05 --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_binary_${trait}
else
       trait=Prostate
       regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --bt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait} &
       wait
       regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE.txt --covarColList age,age2,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_binary_${trait}_pred.list --bt --firth --approx --pThresh 0.05 --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_binary_${trait}
fi