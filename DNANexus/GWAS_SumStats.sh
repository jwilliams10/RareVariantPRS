docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/Clean_Data/all_chr.bed
dx download UKB_PRS:JW/Clean_Data/all_chr.bim
dx download UKB_PRS:JW/Clean_Data/all_chr.fam
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript GWAS_SumStats.R ${1}"