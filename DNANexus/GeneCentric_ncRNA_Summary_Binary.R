rm(list = ls())

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_Analysis/GeneCentric_ncRNA_Summary_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_Analysis/GeneCentric_ncRNA_Summary_Binary.sh -icmd="bash GeneCentric_ncRNA_Summary_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentric_ncRNA/ --priority low --instance-type mem2_ssd1_v2_x2
# done

library(data.table)
trait <- as.numeric(commandArgs(TRUE)[1])

if(trait == 1){
  trait <- "Asthma"
}else if(trait == 2){
  trait <- "T2D"
}else if(trait == 3){
  trait <- "Breast"
}else if(trait == 4){
  trait <- "Prostate"
}else{
  trait <- "CAD"
}

results_ncRNA_genome <- NULL
for (i in 1:22){
  results <- get(load(paste0(trait,"_UKBB_WES_ncRNA_Train",i,".Rdata")))
  file.remove(paste0(trait,"_UKBB_WES_ncRNA_Train",i,".Rdata"))
  results_ncRNA_genome <- c(results_ncRNA_genome, results)
}

results_ncRNA_genome <- do.call("rbind",results_ncRNA_genome)

results_ncRNA_genome <- data.frame(Gene = unlist(results_ncRNA_genome[,c("Gene name")]),
                                   Chr = unlist(results_ncRNA_genome[,c("Chr")]),
                                   Category = unlist(results_ncRNA_genome[,c("Category")]),
                                   Number_SNV = unlist(results_ncRNA_genome[,c("#SNV")]),
                                   Burden_1_1 = unlist(results_ncRNA_genome[,c("Burden(1,1)")]),
                                   STAARO = unlist(results_ncRNA_genome[,c("STAAR-O")]))

results_ncRNA_genome$Number_SNV <- as.numeric(results_ncRNA_genome$Number_SNV)
results_ncRNA_genome <- results_ncRNA_genome[results_ncRNA_genome$Number_SNV < 2000,]

write.csv(results_ncRNA_genome,file = paste0(trait,"_ncRNA_sig.csv"),row.names = FALSE)