# load bigsnpr image
docker load -i r_with_plink.tar.gz

if [ $1 = 1 ]
then
       trait=Asthma
 elif [ $1 = 2 ]
then
       trait=CAD
 elif [ $1 = 3 ]
then
       trait=T2D
 elif [ $1 = 4 ]
then 
       trait=Breast
else
       trait=Prostate
fi


dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/${trait}_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/${trait}_Best_Validation_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}Tune_BestPRS.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}Validation_BestPRS.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/CT/${trait}_prs_validation_best.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_ldpred2_validation_prs_best.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_lassosum2_validation_prs_best.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Common_Plus_Rare_PRS_Binary.R ${1}"