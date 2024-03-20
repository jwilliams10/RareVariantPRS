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

CHR=$(sed "${2}q;d" coding_chr.csv | awk -F ',' '{print $2}')

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/${trait}_coding_sig_es.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Tune_Null_Model.RData
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Validation_Null_Model.RData
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GeneCentric_Coding_PRS.R ${2}"