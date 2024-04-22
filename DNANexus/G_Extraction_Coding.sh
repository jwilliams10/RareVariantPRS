# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/coding_sig.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/noncoding_sig.csv

for trait in BMI HDL LDL Height TC logTG;
do
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Train_Null_Model.RData
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Tune_Null_Model.RData
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Validation_Null_Model.RData
done

dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${1}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript G_Extraction_Coding.R"