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

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/${trait}_UKBB_WES_Noncoding_Train${i}.Rdata
done

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GeneCentric_ncRNA_Summary.R ${1}"