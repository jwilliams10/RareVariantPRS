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

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricNoncoding/${trait}_UKBB_WES_Noncoding_Train${i}.Rdata
done

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GeneCentric_Noncoding_Summary_Binary.R ${1}"