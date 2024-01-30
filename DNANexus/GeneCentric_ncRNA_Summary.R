rm(list = ls())

# for array in 1 2 3 4 5 6;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/GeneCentric_ncRNA_Summary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/GeneCentric_ncRNA_Summary.sh -icmd="bash GeneCentric_ncRNA_Summary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/ --priority low --instance-type mem2_ssd1_v2_x8
# done

library(data.table)
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

results_ncRNA_genome <- NULL
for (i in 1:22){
  results <- get(load(paste0(trait,"_UKBB_WES_Noncoding_Train",i,".Rdata")))
  file.remove(paste0(trait,"_UKBB_WES_Noncoding_Train",i,".Rdata"))
  results_ncRNA_genome <- c(results_ncRNA_genome, results)
}

results_ncRNA_genome <- do.call("rbind",results_ncRNA_genome)

results_ncRNA_genome <- data.frame(Gene = unlist(results_ncRNA_genome[,c("Gene name")]),
                                   Chr = unlist(results_ncRNA_genome[,c("Chr")]),
                                   Category = unlist(results_ncRNA_genome[,c("Category")]),
                                   Number_SNV = unlist(results_ncRNA_genome[,c("#SNV")]),
                                   Burden_1_1 = unlist(results_ncRNA_genome[,c("Burden(1,1)")]),
                                   STAARO = unlist(results_ncRNA_genome[,c("STAAR-O")]))

write.csv(results_ncRNA_genome,file = paste0(trait,"_ncRNA_sig.csv"),row.names = FALSE)