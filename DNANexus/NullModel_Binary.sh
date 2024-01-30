# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Asthma_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Asthma_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Asthma_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/CAD_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/CAD_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/CAD_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/T2D_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/T2D_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/T2D_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Breast_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Breast_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Breast_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Prostate_Best_Train_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Prostate_Best_Tune_All.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/Prostate_Best_Validation_All.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript NullModel_Binary.R ${1}"