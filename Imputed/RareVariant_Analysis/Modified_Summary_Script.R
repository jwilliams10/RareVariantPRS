# =============================================================================
# [RICE-ANNOTATION] Imputed/RareVariant_Analysis/Modified_Summary_Script.R
# Purpose: Analysis/helper script supporting the pipeline in this directory (see README in this folder).
#
# Paper linkage:
#   - Manuscript: Results -> UKB Imputed + WES Results (Fig. 4–5).
#   - Supplementary Figures: Supp. Fig. 3–4 (association diagnostics) and Supp. Fig. 7–11 (PRS performance + sensitivity).
#   - Supplementary Data: key UKB cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list = ls())

library(data.table)

gene_centric_coding_jobs_num <- 381
input_path <- "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/"
output_path <- input_path

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")

binary_traits <- c("Breast","Prostate","CAD","T2D","Asthma")

i <- 1

for(trait in c(continuous_traits,binary_traits)){
  summary_time <- system.time({
    gene_centric_results_name <- paste0(trait,"_UKBB_WES_Coding_Train")
    
    results_coding_genome <- NULL
    for (i in 1:gene_centric_coding_jobs_num){
      results <- get(load(paste0(input_path,gene_centric_results_name,"_",i,".Rdata")))
      results_coding_genome <- c(results_coding_genome, results)
    }
    
    combine_jake <- function(x){
      x <- as.data.frame(x)
      a <- data.frame(Gene = unlist(x$`Gene name`),
                      Chr = unlist(x$Chr),
                      Category = unlist(x$Category),
                      Number_SNV = unlist(x$`#SNV`),
                      Burden_1_1 = unlist(x$`Burden(1,1)`),
                      STAARB = unlist(x$`STAAR-B(1,1)`))
      return(a)
    }
    
    results_coding_genome <- lapply(results_coding_genome, combine_jake)
    results_coding_genome <- rbindlist(results_coding_genome)
    results_coding_genome$Number_SNV <- as.numeric(results_coding_genome$Number_SNV)
    results_coding_genome <- results_coding_genome[results_coding_genome$Number_SNV < 2000,]
    
    write.csv(results_coding_genome,paste0(output_path,trait,"_coding_sig.csv"),row.names = FALSE)
    
    
    analysis_time <- 0
    for (i in 1:gene_centric_coding_jobs_num){
      time <- get(load(paste0(input_path,gene_centric_results_name,"_",i,"_Time.Rdata")))
      analysis_time <- analysis_time + time
    }
    save(analysis_time,file = paste0(output_path,trait,"_Analysis_Time.RData"))
  })[3]
  save(summary_time,file = paste0(output_path,trait,"_Summary_Time.RData"))
}

for(trait in c(continuous_traits,binary_traits)){
  gene_centric_results_name <- paste0(trait,"_UKBB_WES_Coding_Train")
  for (i in 1:gene_centric_coding_jobs_num){
    file.remove(paste0(input_path,gene_centric_results_name,"_",i,".Rdata"))
    file.remove(paste0(input_path,gene_centric_results_name,"_",i,"_Time.Rdata"))
  }
}