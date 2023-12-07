rm(list = ls())

merge_list_train <- paste0("/data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr",2:22,"/ukbb_wes_200k_chr",2:22,"_common")
write.table(merge_list_train,file = "/data/williamsjacr/UKB_WES_Full_Processed_Data/merge_list.txt",col.names = F,row.names = F,quote=F)
system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr",1,"/ukbb_wes_200k_chr",1,"_common --merge-list /data/williamsjacr/UKB_WES_Full_Processed_Data/merge_list.txt --make-bed --out /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr"))

for(i in 1:22){
  system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr",i,"/ukbb_wes_200k_chr",i,"_common --exclude /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr-merge.missnp --make-bed --out /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr",i,"/ukbb_wes_200k_chr",i,"_common"))
}

system(paste0("/data/williamsjacr/software/plink --bfile /data/williamsjacr/UKB_WES_Full_Processed_Data/pVCF/chr",1,"/ukbb_wes_200k_chr",1,"_common --merge-list /data/williamsjacr/UKB_WES_Full_Processed_Data/merge_list.txt --make-bed --out /data/williamsjacr/UKB_WES_Full_Processed_Data/all_chr"))
