rm(list = ls())

# for array in 1 2 3 4 5 6;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/GWAS_SumStats.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/GWAS_SumStats.sh -icmd="bash GWAS_SumStats.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GWAS_SummaryStats/ --instance-type mem1_ssd1_v2_x36
# done

All_Train <- read.delim("All_Train.txt")
All_Train <- subset(All_Train,select = -c(FID))
write.table(All_Train,file = "All_Train.txt",sep = '\t',row.names = FALSE,quote = FALSE)

all_chr <- read.table("all_chr.fam", quote="\"", comment.char="")
all_chr[,1] <- 0  
write.table(all_chr,"all_chr.fam",row.names = FALSE,col.names = FALSE)

trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "BMI"
}else if(trait == 2){
  trait <- "TC"
}else if(trait == 3){
  trait <- "HDL"
}else if(trait == 4){
  trait <- "LDL"
}else if(trait == 5){
  trait <- "logTG"
}else{
  trait <- "Height"
}

system(paste0("plink2 --bfile all_chr --pheno All_Train.txt --pheno-name ",trait," --linear --covar All_Train.txt --covar-name age, age2, sex, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10 --vif 999 --out ",trait,"_sumstats"))

system("rm all_chr.bed")
system("rm all_chr.fam")
system("rm all_chr.bim")
system("rm All_Train.txt")