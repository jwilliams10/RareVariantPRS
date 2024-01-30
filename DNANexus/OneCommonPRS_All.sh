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

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/CT/${trait}_prs_all_train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/CT/${trait}_prs_all_tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/CT/${trait}_prs_all_validation.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_train_ldpred2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_tune_ldpred2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_validation_ldpred2.sscore

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_train_lassosum2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_tune_lassosum2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_validation_lassosum2.sscore

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript OneCommonPRS_All.R ${1}"