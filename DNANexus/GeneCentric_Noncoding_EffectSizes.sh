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

if [ $2 -lt 185 ]
then
       CHR=1
 elif [ $2 -lt 289 ]
then
       CHR=2
 elif [ $2 -lt 385 ]
then
       CHR=3
 elif [ $2 -lt 449 ]
then 
       CHR=4
 elif [ $2 -lt 529 ]
then 
       CHR=5
 elif [ $2 -lt 617 ]
then 
       CHR=6
 elif [ $2 -lt 689 ]
then 
       CHR=7
 elif [ $2 -lt 745 ]
then 
       CHR=8
 elif [ $2 -lt 817 ]
then 
       CHR=9
 elif [ $2 -lt 881 ]
then 
       CHR=10
 elif [ $2 -lt 985 ]
then 
       CHR=11
 elif [ $2 -lt 1073 ]
then 
       CHR=12
 elif [ $2 -lt 1097 ]
then 
       CHR=13
 elif [ $2 -lt 1153 ]
then 
       CHR=14
 elif [ $2 -lt 1201 ]
then 
       CHR=15
 elif [ $2 -lt 1265 ]
then 
       CHR=16
 elif [ $2 -lt 1361 ]
then 
       CHR=17
 elif [ $2 -lt 1385 ]
then 
       CHR=18
 elif [ $2 -lt 1505 ]
then 
       CHR=19
 elif [ $2 -lt 1553 ]
then 
       CHR=20
 elif [ $2 -lt 1569 ]
then 
       CHR=21
else
       CHR=22
fi

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/${trait}_noncoding_sig.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Train_Null_Model.RData
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GeneCentric_Noncoding_EffectSizes.R ${2}"