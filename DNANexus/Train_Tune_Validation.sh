docker load -i r_with_plink.tar.gz

dx download UKB_PRS:JW/Clean_Data/all_chr.fam
dx download UKB_PRS:UKB_200K_WGS_AGDS_uncompressed_newQClabel/ukb.200k.wgs.chr22.pass.annotated.gds
dx download UKB_PRS:ukb_multi.cov
dx download UKB_PRS:ukb_multi.pheno
dx download UKB_PRS:ukb_multi_anc.RDS

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript Train_Tune_Validation.R"