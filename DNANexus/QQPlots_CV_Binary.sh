# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_act_Asthma.regenie
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_act_CAD.regenie
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_act_T2D.regenie
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_bp_Prostate.regenie
dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_bp_Breast.regenie
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript QQPlots_CV_Binary.R ${1}"