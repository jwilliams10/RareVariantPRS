rm(list = ls())

merge_list_train <- paste0("/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",2:22,"_common_train")
merge_list_tune <- paste0("/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",2:22,"_common_tune")
merge_list_validation <- paste0("/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",2:22,"_common_validation")
merge_list_reference <- paste0("/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",2:22,"_common_reference")
write.table(merge_list_train,file = "/data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_train.txt",col.names = F,row.names = F,quote=F)
write.table(merge_list_tune,file = "/data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_tune.txt",col.names = F,row.names = F,quote=F)
write.table(merge_list_validation,file = "/data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_validation.txt",col.names = F,row.names = F,quote=F)
write.table(merge_list_reference,file = "/data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_reference.txt",col.names = F,row.names = F,quote=F)
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",1,"_common_train --merge-list /data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_train.txt --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr"))

for(i in 1:22){
  system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_train --exclude /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr-merge.missnp --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_train"))
  system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_tune --exclude /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr-merge.missnp --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_tune"))
  system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_validation --exclude /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr-merge.missnp --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_validation"))
  system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_reference --exclude /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr-merge.missnp --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",i,"_common_reference"))
}

system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",1,"_common_train --merge-list /data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_train.txt --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_train"))
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",1,"_common_tune --merge-list /data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_tune.txt --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_tune"))
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",1,"_common_validation --merge-list /data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_validation.txt --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_validation"))
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",1,"_common_reference --merge-list /data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list_reference.txt --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr_reference"))
