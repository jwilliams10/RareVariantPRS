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

for trait in BMI TC HDL LDL logTG Height;do

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/${trait}_UKBB_WES_Coding_Train${i}.Rdata
done

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/${trait}_UKBB_WES_Noncoding_Train${i}.Rdata
done

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/${trait}_UKBB_WES_ncRNA_Train${i}.Rdata
done

done

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript RV_Analysis_Summary.R"