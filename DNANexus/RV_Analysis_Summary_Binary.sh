# load bigsnpr image
docker load -i r_with_plink.tar.gz

for trait in Asthma CAD T2D Breast Prostate;do

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricCoding/${trait}_UKBB_WES_Coding_Train${i}.Rdata
done

for i in {1..379};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricNoncoding/${trait}_UKBB_WGS_Noncoding_Train${i}.Rdata
done

for i in {1..22};do
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentric_ncRNA/${trait}_UKBB_WES_ncRNA_Train${i}.Rdata
done

done

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript RV_Analysis_Summary_Binary.R"