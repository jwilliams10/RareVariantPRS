rm(list = ls())


# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/RV_Analysis_Summary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_Analysis/RV_Analysis_Summary.sh -icmd="bash RV_Analysis_Summary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/ --priority low --instance-type mem2_ssd1_v2_x4


library(data.table)

results_coding_genome <- NULL
results_noncoding_genome <- NULL

for(trait in c("BMI","TC","LDL","logTG","Height")){
  results_coding_genome_tmp <- NULL
  for (i in 1:22){
    results <- get(load(paste0(trait,"_UKBB_WES_Coding_Train",i,".Rdata")))
    file.remove(paste0(trait,"_UKBB_WES_Coding_Train",i,".Rdata"))
    results_coding_genome_tmp <- c(results_coding_genome_tmp, results)
  }
  
  combine_jake <- function(x){
    a <- data.frame(Gene = c(unlist(x$plof[,c("Gene name")]),unlist(x$plof_ds[,c("Gene name")]),unlist(x$missense[,c("Gene name")]),unlist(x$disruptive_missense[,c("Gene name")]),unlist(x$synonymous[,c("Gene name")])),
                    Chr = c(unlist(x$plof[,c("Chr")]),unlist(x$plof_ds[,c("Chr")]),unlist(x$missense[,c("Chr")]),unlist(x$disruptive_missense[,c("Chr")]),unlist(x$synonymous[,c("Chr")])),
                    Category = c(unlist(x$plof[,c("Category")]),unlist(x$plof_ds[,c("Category")]),unlist(x$missense[,c("Category")]),unlist(x$disruptive_missense[,c("Category")]),unlist(x$synonymous[,c("Category")])),
                    Number_SNV = c(unlist(x$plof[,c("#SNV")]),unlist(x$plof_ds[,c("#SNV")]),unlist(x$missense[,c("#SNV")]),unlist(x$disruptive_missense[,c("#SNV")]),unlist(x$synonymous[,c("#SNV")])),
                    Burden_1_1 = c(unlist(x$plof[,c("Burden(1,1)")]),unlist(x$plof_ds[,c("Burden(1,1)")]),unlist(x$missense[,c("Burden(1,1)")]),unlist(x$disruptive_missense[,c("Burden(1,1)")]),unlist(x$synonymous[,c("Burden(1,1)")])),
                    STAARB = c(unlist(x$plof[,c("STAAR-B(1,1)")]),unlist(x$plof_ds[,c("STAAR-B(1,1)")]),unlist(x$missense[,c("STAAR-B(1,1)")]),unlist(x$disruptive_missense[,c("STAAR-B(1,1)")]),unlist(x$synonymous[,c("STAAR-B(1,1)")])))
    return(a)
  }
  
  results_coding_genome_tmp <- lapply(results_coding_genome_tmp, combine_jake)
  results_coding_genome_tmp <- rbindlist(results_coding_genome_tmp)
  results_coding_genome_tmp$Number_SNV <- as.numeric(results_coding_genome_tmp$Number_SNV)
  results_coding_genome_tmp <- results_coding_genome_tmp[results_coding_genome_tmp$Number_SNV < 2000,]
  
  results_coding_genome_tmp$Trait <- trait
  
  results_coding_genome <- rbind(results_coding_genome,results_coding_genome_tmp)
  
  
  results_noncoding_genome_tmp <- NULL
  for (i in 1:22){
    results <- get(load(paste0(trait,"_UKBB_WES_Noncoding_Train",i,".Rdata")))
    file.remove(paste0(trait,"_UKBB_WES_Noncoding_Train",i,".Rdata"))
    results_noncoding_genome_tmp <- c(results_noncoding_genome_tmp, results)
  }
  
  combine_jake <- function(x){
    a <- data.frame(Gene = c(unlist(x$upstream[,c("Gene name")]),unlist(x$downstream[,c("Gene name")]),unlist(x$UTR[,c("Gene name")]),unlist(x$promoter_CAGE[,c("Gene name")]),unlist(x$promoter_DHS[,c("Gene name")]),unlist(x$enhancer_CAGE[,c("Gene name")]),unlist(x$enhancer_DHS[,c("Gene name")])),
                    Chr = c(unlist(x$upstream[,c("Chr")]),unlist(x$downstream[,c("Chr")]),unlist(x$UTR[,c("Chr")]),unlist(x$promoter_CAGE[,c("Chr")]),unlist(x$promoter_DHS[,c("Chr")]),unlist(x$enhancer_CAGE[,c("Chr")]),unlist(x$enhancer_DHS[,c("Chr")])),
                    Category = c(unlist(x$upstream[,c("Category")]),unlist(x$downstream[,c("Category")]),unlist(x$UTR[,c("Category")]),unlist(x$promoter_CAGE[,c("Category")]),unlist(x$promoter_DHS[,c("Category")]),unlist(x$enhancer_CAGE[,c("Category")]),unlist(x$enhancer_DHS[,c("Category")])),
                    Number_SNV = c(unlist(x$upstream[,c("#SNV")]),unlist(x$downstream[,c("#SNV")]),unlist(x$UTR[,c("#SNV")]),unlist(x$promoter_CAGE[,c("#SNV")]),unlist(x$promoter_DHS[,c("#SNV")]),unlist(x$enhancer_CAGE[,c("#SNV")]),unlist(x$enhancer_DHS[,c("#SNV")])),
                    Burden_1_1 = c(unlist(x$upstream[,c("Burden(1,1)")]),unlist(x$downstream[,c("Burden(1,1)")]),unlist(x$UTR[,c("Burden(1,1)")]),unlist(x$promoter_CAGE[,c("Burden(1,1)")]),unlist(x$promoter_DHS[,c("Burden(1,1)")]),unlist(x$enhancer_CAGE[,c("Burden(1,1)")]),unlist(x$enhancer_DHS[,c("Burden(1,1)")])),
                    STAARB = c(unlist(x$upstream[,c("STAAR-B(1,1)")]),unlist(x$downstream[,c("STAAR-B(1,1)")]),unlist(x$UTR[,c("STAAR-B(1,1)")]),unlist(x$promoter_CAGE[,c("STAAR-B(1,1)")]),unlist(x$promoter_DHS[,c("STAAR-B(1,1)")]),unlist(x$enhancer_CAGE[,c("STAAR-B(1,1)")]),unlist(x$enhancer_DHS[,c("STAAR-B(1,1)")])))
    return(a)
  }
  
  results_noncoding_genome_tmp <- lapply(results_noncoding_genome_tmp, combine_jake)
  results_noncoding_genome_tmp <- rbindlist(results_noncoding_genome_tmp)
  
  results_noncoding_genome_tmp$Number_SNV <- as.numeric(results_noncoding_genome_tmp$Number_SNV)
  results_noncoding_genome_tmp <- results_noncoding_genome_tmp[results_noncoding_genome_tmp$Number_SNV < 2000,]
  
  
  
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
                                     STAARB = unlist(results_ncRNA_genome[,c("STAAR-B(1,1)")]))
  
  results_ncRNA_genome$Number_SNV <- as.numeric(results_ncRNA_genome$Number_SNV)
  results_ncRNA_genome <- results_ncRNA_genome[results_ncRNA_genome$Number_SNV < 2000,]
  
  results_noncoding_genome_tmp <- cbind(results_noncoding_genome_tmp,results_ncRNA_genome)
  results_noncoding_genome_tmp$Trait <- trait
  
  results_noncoding_genome <- rbind(results_noncoding_genome,results_noncoding_genome_tmp)
  
}

results_coding_genome <- results_coding_genome[results_coding_genome$STAARB < 1e-3,]
results_noncoding_genome <- results_noncoding_genome[results_noncoding_genome$STAARB < 1e-3,]

write.csv(results_coding_genome,file = paste0("coding_sig.csv"),row.names = FALSE)
write.csv(results_noncoding_genome,file = paste0("noncoding_sig.csv"),row.names = FALSE)

system("rm *.Rdata")