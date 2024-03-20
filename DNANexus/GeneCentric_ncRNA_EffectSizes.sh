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

if [ $2 -lt 19 ]
then
       CHR=1
 elif [ $2 -lt 41 ]
then
       CHR=2
 elif [ $2 -lt 56 ]
then
       CHR=3
 elif [ $2 -lt 65 ]
then 
       CHR=4
 elif [ $2 -lt 76 ]
then 
       CHR=5
 elif [ $2 -lt 85 ]
then 
       CHR=6
 elif [ $2 -lt 94 ]
then 
       CHR=7
 elif [ $2 -lt 102 ]
then 
       CHR=8
 elif [ $2 -lt 110 ]
then 
       CHR=9
 elif [ $2 -lt 117 ]
then 
       CHR=10
 elif [ $2 -lt 126 ]
then 
       CHR=11
 elif [ $2 -lt 134 ]
then 
       CHR=12
 elif [ $2 -lt 143 ]
then 
       CHR=13
 elif [ $2 -lt 150 ]
then 
       CHR=14
 elif [ $2 -lt 158 ]
then 
       CHR=15
 elif [ $2 -lt 165 ]
then 
       CHR=16
 elif [ $2 -lt 174 ]
then 
       CHR=17
 elif [ $2 -lt 179 ]
then 
       CHR=18
 elif [ $2 -lt 187 ]
then 
       CHR=19
 elif [ $2 -lt 193 ]
then 
       CHR=20
 elif [ $2 -lt 198 ]
then 
       CHR=21
else
       CHR=22
fi

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/${trait}_ncRNA_sig.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Train_Null_Model.RData
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GeneCentric_ncRNA_EffectSizes.R ${2}"