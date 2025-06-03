#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --array=1-6
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5G

if [ ${SLURM_ARRAY_TASK_ID} = 1 ]
then
       trait=BMI
 elif [ ${SLURM_ARRAY_TASK_ID} = 2 ]
then
       trait=TC
 elif [ ${SLURM_ARRAY_TASK_ID} = 3 ]
then
       trait=HDL
 elif [ ${SLURM_ARRAY_TASK_ID} = 4 ]
then 
       trait=LDL
 elif [ ${SLURM_ARRAY_TASK_ID} = 5 ]
then 
       trait=logTG
else
       trait=Height
fi

# module purge
module load regenie/3.0.3

regenie --step 1 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega_regenie_step1 --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE_NewPCs.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE_NewPCs.txt --covarColList age,age2,sex,NewPCs1,NewPCs2,NewPCs3,NewPCs4,NewPCs5,NewPCs6,NewPCs7,NewPCs8,NewPCs9,NewPCs10 --qt --bsize 1000 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_continuous_NewPCs_${trait} &
wait
regenie --step 2 --bed /data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/ukb_hm3_mega --phenoFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE_NewPCs.txt --phenoColList ${trait} --covarFile /data/williamsjacr/UKB_WES_Phenotypes/All_Train_${trait}_REGENIE_NewPCs.txt --covarColList age,age2,sex,NewPCs1,NewPCs2,NewPCs3,NewPCs4,NewPCs5,NewPCs6,NewPCs7,NewPCs8,NewPCs9,NewPCs10 --pred /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step1_continuous_NewPCs_${trait}_pred.list --qt --bsize 400 --out /data/williamsjacr/UKB_WES_Phenotypes/Imputed/GWAS_Summary_Statistics/regenie_step2_continuous_NewPCs