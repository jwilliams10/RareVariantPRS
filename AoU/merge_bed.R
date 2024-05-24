# %%writefile score_task.R
# 
# tasks <- data.frame(check.names = FALSE)
# 
# tasks <- rbind(tasks, data.frame(
#   '--input R_Script'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Scripts/merge_bed.R",
#   '--input-recursive INPUT_PATH'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants",
#   '--output-recursive OUTPUT_PATH'="gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Data/CommonVariants",
#   # '--output out_file'=paste0("gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/JW/AoU_Phenotypes/Results/test",chrom,".Rdata"),
#   check.names = FALSE
# ))
# 
# write.table(tasks, 
#             file="score_task.txt", 
#             row.names=F, col.names=T, 
#             sep='\t', quote=F)


%%writefile merge_bed.R

rm(list = ls())

INPUT_PATH <- commandArgs(TRUE)[1]
print(INPUT_PATH)

OUTPUT_PATH <- commandArgs(TRUE)[2]
print(OUTPUT_PATH)

merge_list_train <- paste0(INPUT_PATH,"/chr",2:22,"_corrected")
write.table(merge_list_train,file = "merge_list.txt",col.names = F,row.names = F,quote=F)
system(paste0("plink --bfile ",INPUT_PATH,"/chr",1,"_corrected --merge-list merge_list.txt --make-bed --out ",OUTPUT_PATH,"/all_chr"))