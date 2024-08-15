# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download -r UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript QQ_Plots_CV_Continuous.R ${1}"