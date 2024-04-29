# load bigsnpr image
docker load -i r_with_plink.tar.gz

for trait in BMI TC HDL LDL logTG Height;do

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/${trait}_UKBB_WES_Coding_Train${i}.Rdata
done

for i in {1..379};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/${trait}_UKBB_WGS_Noncoding_Train_${i}.Rdata
done

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/${trait}_UKBB_WES_ncRNA_Train${i}.Rdata
done

done

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript RV_Analysis_Summary.R"