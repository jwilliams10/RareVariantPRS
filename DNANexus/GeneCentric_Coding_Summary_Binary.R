rm(list = ls())

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_Analysis/GeneCentric_Coding_Summary_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_Analysis/GeneCentric_Coding_Summary_Binary.sh -icmd="bash GeneCentric_Coding_Summary_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricCoding/ --priority low --instance-type mem2_ssd1_v2_x2
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

results_coding_genome <- NULL
for (i in 1:22){
  results <- get(load(paste0(trait,"_UKBB_WES_Coding_Train",i,".Rdata")))
  file.remove(paste0(trait,"_UKBB_WES_Coding_Train",i,".Rdata"))
  results_coding_genome <- c(results_coding_genome, results)
}

combine_jake <- function(x){
  a <- data.frame(Gene = c(unlist(x$plof[,c("Gene name")]),unlist(x$plof_ds[,c("Gene name")]),unlist(x$missense[,c("Gene name")]),unlist(x$disruptive_missense[,c("Gene name")]),unlist(x$synonymous[,c("Gene name")])),
                  Chr = c(unlist(x$plof[,c("Chr")]),unlist(x$plof_ds[,c("Chr")]),unlist(x$missense[,c("Chr")]),unlist(x$disruptive_missense[,c("Chr")]),unlist(x$synonymous[,c("Chr")])),
                  Category = c(unlist(x$plof[,c("Category")]),unlist(x$plof_ds[,c("Category")]),unlist(x$missense[,c("Category")]),unlist(x$disruptive_missense[,c("Category")]),unlist(x$synonymous[,c("Category")])),
                  Number_SNV = c(unlist(x$plof[,c("#SNV")]),unlist(x$plof_ds[,c("#SNV")]),unlist(x$missense[,c("#SNV")]),unlist(x$disruptive_missense[,c("#SNV")]),unlist(x$synonymous[,c("#SNV")])),
                  Burden_1_1 = c(unlist(x$plof[,c("Burden(1,1)")]),unlist(x$plof_ds[,c("Burden(1,1)")]),unlist(x$missense[,c("Burden(1,1)")]),unlist(x$disruptive_missense[,c("Burden(1,1)")]),unlist(x$synonymous[,c("Burden(1,1)")])),
                  STAARO = c(unlist(x$plof[,c("STAAR-O")]),unlist(x$plof_ds[,c("STAAR-O")]),unlist(x$missense[,c("STAAR-O")]),unlist(x$disruptive_missense[,c("STAAR-O")]),unlist(x$synonymous[,c("STAAR-O")])))
  return(a)
}

results_coding_genome <- lapply(results_coding_genome, combine_jake)
results_coding_genome <- rbindlist(results_coding_genome)
results_coding_genome$Number_SNV <- as.numeric(results_coding_genome$Number_SNV)
results_coding_genome <- results_coding_genome[results_coding_genome$Number_SNV < 2000,]

write.csv(results_coding_genome,file = paste0(trait,"_coding_sig.csv"),row.names = FALSE)