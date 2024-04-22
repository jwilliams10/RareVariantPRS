# load bigsnpr image
docker load -i r_with_plink.tar.gz

for trait in BMI HDL LDL Height TC logTG;
do
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}_coding_sig.csv
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}_noncoding_sig.csv
  
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Train_Null_Model.RData
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Tune_Null_Model.RData
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Validation_Null_Model.RData
done

for CHR in {1..22};
do
  dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
done 

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Single_RareVariant_PRS_All.R"