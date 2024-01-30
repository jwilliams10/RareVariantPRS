rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/merge_bed.R -icmd="Rscript merge_bed.R" --destination UKB_PRS:JW/ --instance-type mem1_ssd1_v2_x36

system("dx ls")
system("dx download UKB_PRS:/JW/Clean_Data/ -r")

merge_list_train <- paste0("Clean_Data/chr",2:22,"_filtered_common")
write.table(merge_list_train,file = "Clean_Data/merge_list.txt",col.names = F,row.names = F,quote=F)
system(paste0("plink --bfile Clean_Data/chr",1,"_filtered_common --merge-list Clean_Data/merge_list.txt --make-bed --out Clean_Data/all_chr"))

# for(i in 1:22){
#   system(paste0("plink --bfile Clean_Data/chr",i,"_filtered_common --exclude Clean_Data/all_chr-merge.missnp --make-bed --out Clean_Data/chr",i,"_filtered_common"))
#   system(paste0("rm Clean_Data/chr",i,"_filtered_common.bed~"))
#   system(paste0("rm Clean_Data/chr",i,"_filtered_common.fam~"))
#   system(paste0("rm Clean_Data/chr",i,"_filtered_common.bim~"))
# }
# 
# system(paste0("plink --bfile Clean_Data/chr",1,"_filtered_common --merge-list Clean_Data/merge_list.txt --make-bed --out Clean_Data/all_chr"))

# system("dx upload /Clean_Data/ -r --path UKB_PRS:JW/")
# 
# system("rm -r Clean_Data")