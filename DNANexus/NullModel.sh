# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/BMI_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/BMI_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/BMI_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/TC_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/TC_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/TC_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/LDL_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/LDL_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/LDL_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/HDL_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/HDL_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/HDL_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/logTG_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/logTG_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/logTG_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/Height_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/Height_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/Height_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript NullModel.R ${1}"