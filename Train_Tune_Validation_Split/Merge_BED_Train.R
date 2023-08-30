rm(list = ls())

merge_list <- paste0("/data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",2:22,"_common_train")
write.table(merge_list,file = "/data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list.txt",col.names = F,row.names = F,quote=F)
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_lipids/Data/split_bed/ukbb_wes_200k_chr",1,"_common_train --merge-list /data/williamsjacr/UKB_WES_lipids/Data/split_bed/merge_list.txt --make-bed --out /data/williamsjacr/UKB_WES_lipids/Data/split_bed/all_chr"))