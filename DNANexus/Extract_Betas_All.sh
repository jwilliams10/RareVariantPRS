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


dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/CT/${trait}_Final_Coefficients.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/CT/${trait}_prs_all_validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/CT/${trait}_prs_all_tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_Final_Coefficients_LASSOSum.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_validation_lassosum2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_tune_lassosum2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_Final_Coefficients_LDPred2.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_validation_ldpred2.sscore
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}_prs_tune_ldpred2.sscore

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/${trait}_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/${trait}_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/validation.txt
dx download UKB_PRS:JW/Clean_Data/all_chr.bed
dx download UKB_PRS:JW/Clean_Data/all_chr.bim
dx download UKB_PRS:JW/Clean_Data/all_chr.fam


# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Extract_Betas_All.R ${1}"