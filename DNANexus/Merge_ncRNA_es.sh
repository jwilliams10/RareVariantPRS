# load bigsnpr image
docker load -i r_with_plink.tar.gz

dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA_es/ -r

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Merge_ncRNA_es.R"