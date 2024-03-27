# load bigsnpr image
docker load -i r_with_plink.tar.gz

for trait in Asthma CAD T2D Breast Prostate;
do
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/${trait}_Best_Tune_All.txt
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/${trait}_Best_Validation_All.txt
  
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}_Best_All_STAARO_Tune_All.txt
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}_Best_All_STAARO_Validation_All.txt
  
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}_Best_All_Burden_Tune_All.txt
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}_Best_All_Burden_Validation_All.txt
done

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Common_Plus_Rare_PRS_Binary.R"