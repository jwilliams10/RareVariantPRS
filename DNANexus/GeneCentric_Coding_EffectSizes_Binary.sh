# load bigsnpr image
docker load -i r_with_plink.tar.gz

if [ $1 = 1 ]
then
       trait=Asthma
 elif [ $1 = 2 ]
then
       trait=T2D
 elif [ $1 = 3 ]
then
       trait=Breast
 elif [ $1 = 4 ]
then 
       trait=Prostate
else
       trait=CAD
fi

if [ $2 -lt 45 ]
then
       CHR=1
 elif [ $2 -lt 71 ]
then
       CHR=2
 elif [ $2 -lt 95 ]
then
       CHR=3
 elif [ $2 -lt 111 ]
then 
       CHR=4
 elif [ $2 -lt 129 ]
then 
       CHR=5
 elif [ $2 -lt 151 ]
then 
       CHR=6
 elif [ $2 -lt 171 ]
then 
       CHR=7
 elif [ $2 -lt 185 ]
then 
       CHR=8
 elif [ $2 -lt 201 ]
then 
       CHR=9
 elif [ $2 -lt 217 ]
then 
       CHR=10
 elif [ $2 -lt 245 ]
then 
       CHR=11
 elif [ $2 -lt 267 ]
then 
       CHR=12
 elif [ $2 -lt 273 ]
then 
       CHR=13
 elif [ $2 -lt 287 ]
then 
       CHR=14
 elif [ $2 -lt 299 ]
then 
       CHR=15
 elif [ $2 -lt 317 ]
then 
       CHR=16
 elif [ $2 -lt 341 ]
then 
       CHR=17
 elif [ $2 -lt 347 ]
then 
       CHR=18
 elif [ $2 -lt 379 ]
then 
       CHR=19
 elif [ $2 -lt 389 ]
then 
       CHR=20
 elif [ $2 -lt 393 ]
then 
       CHR=21
else
       CHR=22
fi

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricCoding/${trait}_coding_sig.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/nullmodels_staar/${trait}_Train_Null_Model.RData
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr${CHR}.pass.annotated.gds
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/Annotation_name_catalog.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GeneCentric_Coding_EffectSizes_Binary.R ${2}"