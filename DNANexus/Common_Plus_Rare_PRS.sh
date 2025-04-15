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

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/${trait}_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/${trait}_Best_Validation_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}Tune_BestPRS.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}Validation_BestPRS.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Common_Plus_Rare_PRS.R ${1}"