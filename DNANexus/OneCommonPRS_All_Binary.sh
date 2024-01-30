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

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/CT/${trait}_prs_all_train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/CT/${trait}_prs_all_tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/CT/${trait}_prs_all_validation.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_prs_train_ldpred2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_prs_tune_ldpred2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_prs_validation_ldpred2.sscore

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_prs_train_lassosum2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_prs_tune_lassosum2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/${trait}_prs_validation_lassosum2.sscore

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript OneCommonPRS_All_Binary.R ${1}"