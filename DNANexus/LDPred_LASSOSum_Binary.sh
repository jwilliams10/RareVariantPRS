docker load -i r_with_plink.tar.gz

if [ $1 = 1 ]
then
       trait=Asthma
       dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_act_Asthma.regenie
 elif [ $1 = 2 ]
then
       trait=CAD
       dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_act_CAD.regenie
 elif [ $1 = 3 ]
then
       trait=T2D
       dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_act_T2D.regenie
 elif [ $1 = 4 ]
then 
       trait=Breast
       dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_bp_Breast.regenie
else
       trait=Prostate
       dx download UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GWAS_SummaryStats/regenie_step2_bp_Prostate.regenie
fi

dx download UKB_PRS:JW/Clean_Data/w_hm3.snplist
dx download UKB_PRS:JW/Clean_Data/SNP_GRCh37_38_match_update.rds
dx download UKB_PRS:JW/UKB_Phenotypes/Data/all_chr_reference.rds
dx download UKB_PRS:JW/UKB_Phenotypes/Data/all_chr_reference.bk
dx download UKB_PRS:JW/UKB_Phenotypes/Data/LDPred2_Genetic_Mappings.zip

dx download UKB_PRS:JW/Clean_Data/all_chr.bed
dx download UKB_PRS:JW/Clean_Data/all_chr.bim
dx download UKB_PRS:JW/Clean_Data/all_chr.fam

dx download UKB_PRS:JW/UKB_Phenotypes/Results/train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/validation.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Train.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Tune.txt
dx download UKB_PRS:JW/UKB_Phenotypes/Results/All_Validation.txt

dx download UKB_PRS:JW/UKB_Phenotypes/Results/all_phenotypes.RData

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript LDPred_LASSOSum_Binary.R ${1}"