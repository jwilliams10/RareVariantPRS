rm(list = ls())

# for array in 1 2 3 4 5;
# do
# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_Analysis/GeneCentric_Noncoding_Summary_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_Analysis/GeneCentric_Noncoding_Summary_Binary.sh -icmd="bash GeneCentric_Noncoding_Summary_Binary.sh ${array}" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricNoncoding/ --priority low --instance-type mem2_ssd1_v2_x8
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

results_noncoding_genome <- NULL
for (i in 1:22){
  results <- get(load(paste0(trait,"_UKBB_WES_Noncoding_Train",i,".Rdata")))
  file.remove(paste0(trait,"_UKBB_WES_Noncoding_Train",i,".Rdata"))
  results_noncoding_genome <- c(results_noncoding_genome, results)
}

combine_jake <- function(x){
  a <- data.frame(Gene = c(unlist(x$upstream[,c("Gene name")]),unlist(x$downstream[,c("Gene name")]),unlist(x$UTR[,c("Gene name")]),unlist(x$promoter_CAGE[,c("Gene name")]),unlist(x$promoter_DHS[,c("Gene name")]),unlist(x$enhancer_CAGE[,c("Gene name")]),unlist(x$enhancer_DHS[,c("Gene name")])),
                  Chr = c(unlist(x$upstream[,c("Chr")]),unlist(x$downstream[,c("Chr")]),unlist(x$UTR[,c("Chr")]),unlist(x$promoter_CAGE[,c("Chr")]),unlist(x$promoter_DHS[,c("Chr")]),unlist(x$enhancer_CAGE[,c("Chr")]),unlist(x$enhancer_DHS[,c("Chr")])),
                  Category = c(unlist(x$upstream[,c("Category")]),unlist(x$downstream[,c("Category")]),unlist(x$UTR[,c("Category")]),unlist(x$promoter_CAGE[,c("Category")]),unlist(x$promoter_DHS[,c("Category")]),unlist(x$enhancer_CAGE[,c("Category")]),unlist(x$enhancer_DHS[,c("Category")])),
                  Number_SNV = c(unlist(x$upstream[,c("#SNV")]),unlist(x$downstream[,c("#SNV")]),unlist(x$UTR[,c("#SNV")]),unlist(x$promoter_CAGE[,c("#SNV")]),unlist(x$promoter_DHS[,c("#SNV")]),unlist(x$enhancer_CAGE[,c("#SNV")]),unlist(x$enhancer_DHS[,c("#SNV")])),
                  Burden_1_1 = c(unlist(x$upstream[,c("Burden(1,1)")]),unlist(x$downstream[,c("Burden(1,1)")]),unlist(x$UTR[,c("Burden(1,1)")]),unlist(x$promoter_CAGE[,c("Burden(1,1)")]),unlist(x$promoter_DHS[,c("Burden(1,1)")]),unlist(x$enhancer_CAGE[,c("Burden(1,1)")]),unlist(x$enhancer_DHS[,c("Burden(1,1)")])),
                  STAARO = c(unlist(x$upstream[,c("STAAR-O")]),unlist(x$downstream[,c("STAAR-O")]),unlist(x$UTR[,c("STAAR-O")]),unlist(x$promoter_CAGE[,c("STAAR-O")]),unlist(x$promoter_DHS[,c("STAAR-O")]),unlist(x$enhancer_CAGE[,c("STAAR-O")]),unlist(x$enhancer_DHS[,c("STAAR-O")])))
  return(a)
}

results_noncoding_genome <- lapply(results_noncoding_genome, combine_jake)
results_noncoding_genome <- rbindlist(results_noncoding_genome)

results_noncoding_genome$Number_SNV <- as.numeric(results_noncoding_genome$Number_SNV)
results_noncoding_genome <- results_noncoding_genome[results_noncoding_genome$Number_SNV < 2000,]

write.csv(results_noncoding_genome,file = paste0(trait,"_noncoding_sig.csv"),row.names = FALSE)