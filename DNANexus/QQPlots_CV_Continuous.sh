# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/BMI_sumstats.BMI.glm.linear
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/LDL_sumstats.LDL.glm.linear
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/HDL_sumstats.HDL.glm.linear
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/TC_sumstats.TC.glm.linear
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/logTG_sumstats.logTG.glm.linear
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/Height_sumstats.Height.glm.linear
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript QQPlots_CV_Continuous.R ${1}"