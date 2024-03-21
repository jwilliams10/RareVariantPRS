# load bigsnpr image
docker load -i r_with_plink.tar.gz

for trait in BMI HDL LDL Height TC logTG;
do
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/${trait}_Noncoding_Burden_PRS.csv
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/${trait}_Noncoding_STAARO_PRS.csv

  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/${trait}_Coding_Burden_PRS.csv
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/${trait}_Coding_STAARO_PRS.csv
  
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/${trait}_ncRNA_Burden_PRS.csv
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/${trait}_ncRNA_STAARO_PRS.csv
  
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Tune_Null_Model.RData
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Validation_Null_Model.RData
done

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Single_RareVariant_PRS_All.R"