# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/coding_sig.csv
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/noncoding_sig.csv

if [ $1 = 1 ]
then
       trait=Asthma
 elif [ $1 = 2 ]
then
       trait=CAD
 elif [ $1 = 3 ]
then
       trait=T2D
 elif [ $1 = 4 ]
then 
       trait=Breast
else
       trait=Prostate
fi

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/nullmodels_staar/${trait}_Train_Null_Model.RData
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/nullmodels_staar/${trait}_Tune_Null_Model.RData
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/nullmodels_staar/${trait}_Validation_Null_Model.RData

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/${trait}_BestPRS.csv

dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

for i in {1..22};do
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/G_star/${trait}_G_Star_Coding_Chr${i}.csv
  dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/G_star/${trait}_G_Star_Noncoding_Chr${i}.csv
done


# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript RareVariantCoefficients_Binary.R ${1}"