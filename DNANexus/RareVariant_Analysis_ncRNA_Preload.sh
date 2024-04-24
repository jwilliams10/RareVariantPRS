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
 elif [ $2 -lt 35 ]
then
       CHR=2
 elif [ $2 -lt 46 ]
then
       CHR=3
 elif [ $2 -lt 56 ]
then 
       CHR=4
 elif [ $2 -lt 69 ]
then 
       CHR=5
 elif [ $2 -lt 79 ]
then 
       CHR=6
 elif [ $2 -lt 89 ]
then 
       CHR=7
 elif [ $2 -lt 100 ]
then 
       CHR=8
 elif [ $2 -lt 108 ]
then 
       CHR=9
 elif [ $2 -lt 117 ]
then 
       CHR=10
 elif [ $2 -lt 128 ]
then 
       CHR=11
 elif [ $2 -lt 141 ]
then 
       CHR=12
 elif [ $2 -lt 147 ]
then 
       CHR=13
 elif [ $2 -lt 156 ]
then 
       CHR=14
 elif [ $2 -lt 166 ]
then 
       CHR=15
 elif [ $2 -lt 178 ]
then 
       CHR=16
 elif [ $2 -lt 191 ]
then 
       CHR=17
 elif [ $2 -lt 198 ]
then 
       CHR=18
 elif [ $2 -lt 208 ]
then 
       CHR=19
 elif [ $2 -lt 214 ]
then 
       CHR=20
 elif [ $2 -lt 218 ]
then 
       CHR=21
else
       CHR=22
fi

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/nullmodels_staar/${trait}_Train_Null_Model.RData
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript RareVariant_Analysis_ncRNA_Preload.R ${1} ${2}"