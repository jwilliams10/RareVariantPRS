# load bigsnpr image
docker load -i r_with_plink.tar.gz

if [ $1 = 1 ]
then
       trait=BMI
 elif [ $1 = 2 ]
then
       trait=TC
 elif [ $1 = 3 ]
then
       trait=HDL
 elif [ $1 = 4 ]
then 
       trait=LDL
 elif [ $1 = 5 ]
then 
       trait=logTG
else
       trait=Height
fi

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/${trait}_sumstats.${trait}.glm.linear
dx download UKB_PRS:JW/Clean_Data/w_hm3.snplist
dx download UKB_PRS:JW/Clean_Data/SNP_GRCh37_38_match_update.rds
dx download UKB_PRS:JW/UKB_Phenotypes/Data/all_chr_reference.rds
dx download UKB_PRS:JW/UKB_Phenotypes/Data/all_chr_reference.bk
dx download UKB_PRS:JW/UKB_Phenotypes/Data/LDPred2_Genetic_Mappings.zip

dx download UKB_PRS:JW/Clean_Data/all_chr.bed
dx download UKB_PRS:JW/Clean_Data/all_chr.bim
dx download UKB_PRS:JW/Clean_Data/all_chr.fam

dx download UKB_PRS:JW/UKB_Phenotypes/Results/train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/validation.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript LDPred_LASSOSum.R ${1}"