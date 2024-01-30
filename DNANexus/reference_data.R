rm(list = ls())
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/reference_data.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/reference_data.sh  -icmd="bash reference_data.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Data/ --instance-type mem1_ssd1_v2_x36

if(!("remotes" %in% rownames(installed.packages()))){
  install.packages("remotes",quiet = TRUE)
}

if(!("bigsnpr" %in% rownames(installed.packages()))){
  install.packages("bigsnpr",quiet = TRUE)
}

library(bigsnpr)

for(i in 1:22){
  system(paste0("plink2 --bfile Clean_Data/chr",i,"_filtered_common --keep reference.txt --make-bed --out chr",i,"_filtered_common_reference"))
}

system(paste0("plink2 --bfile Clean_Data/all_chr --keep reference.txt --make-bed --out all_chr_reference"))

bigsnpr::snp_readBed("all_chr_reference.bed",backingfile = "all_chr_reference")

system("rm -r Clean_Data/")
file.remove("reference.txt")