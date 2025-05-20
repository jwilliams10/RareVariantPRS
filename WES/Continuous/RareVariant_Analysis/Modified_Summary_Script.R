rm(list = ls())

library(data.table)

gene_centric_coding_jobs_num <- 381
input_path <- "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/"
output_path <- input_path

trait <- "BMI"
i <- 1

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  gene_centric_results_name <- paste0(trait,"_UKBB_WES_Coding_Train")
  
  results_coding_genome <- NULL
  for (i in 1:gene_centric_coding_jobs_num){
    results <- get(load(paste0(input_path,gene_centric_results_name,"_",i,".Rdata")))
    # system(paste0("rm ",paste0(input_path,gene_centric_results_name,"_",i,".Rdata")))
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
  
  write.csv(results_coding_genome,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_coding_sig.csv"),row.names = FALSE)
}

system(paste0("rm ",paste0(input_path,"*.Rdata")))




# rm(list = ls())
# 
# ncRNA_jobs_num <- 223
# ncRNA_input_path <- "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentric_ncRNA/"
# 
# input_path <- "/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/"
# output_path <- input_path
# ## number of jobs
# gene_centric_noncoding_jobs_num <- 387
# 
# trait <- "BMI"
# i <- 1
# 
# for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
#   gene_centric_results_name <- paste0(trait,"_UKBB_WES_Noncoding_Train")
#   ncRNA_results_name <- paste0(trait,"_UKBB_WES_ncRNA_Train")
#   
#   results_ncRNA_genome <- NULL
#   for (i in 1:ncRNA_jobs_num){
#     results <- get(load(paste0(ncRNA_input_path,ncRNA_results_name,"_",i,".Rdata")))
#     results <- as.data.frame(results)
#     results_ncRNA_genome <- rbind(results_ncRNA_genome, results)
#   }
#   
#   results_ncRNA_genome <- data.frame(Gene = unlist(results_ncRNA_genome[,c("Gene name")]),
#                                      Chr = unlist(results_ncRNA_genome[,c("Chr")]),
#                                      Category = unlist(results_ncRNA_genome[,c("Category")]),
#                                      Number_SNV = unlist(results_ncRNA_genome[,c("#SNV")]),
#                                      Burden_1_1 = unlist(results_ncRNA_genome[,c("Burden(1,1)")]),
#                                      STAARB = unlist(results_ncRNA_genome[,c("STAAR-B(1,1)")]))
#   
#   results_ncRNA_genome$Number_SNV <- as.numeric(results_ncRNA_genome$Number_SNV)
#   results_ncRNA_genome <- results_ncRNA_genome[results_ncRNA_genome$Number_SNV < 2000,]
#   
#   
#   
#   
#   results_noncoding_genome <- NULL
#   for (i in 1:gene_centric_noncoding_jobs_num){
#     results <- get(load(paste0(input_path,gene_centric_results_name,"_",i,".Rdata")))
#     results_noncoding_genome <- c(results_noncoding_genome, results)
#   }
#   
#   combine_jake <- function(x){
#     x <- as.data.frame(x)
#     a <- data.frame(Gene = unlist(x$`Gene name`),
#                     Chr = unlist(x$Chr),
#                     Category = unlist(x$Category),
#                     Number_SNV = unlist(x$`#SNV`),
#                     Burden_1_1 = unlist(x$`Burden(1,1)`),
#                     STAARB = unlist(x$`STAAR-B(1,1)`))
#     return(a)
#   }
#   
#   results_noncoding_genome <- lapply(results_noncoding_genome, combine_jake)
#   results_noncoding_genome <- rbindlist(results_noncoding_genome)
#   
#   results_noncoding_genome$Number_SNV <- as.numeric(results_noncoding_genome$Number_SNV)
#   results_noncoding_genome <- results_noncoding_genome[results_noncoding_genome$Number_SNV < 2000,]
#   
#   
#   
#   
#   
#   results_noncoding_genome <- rbind(results_noncoding_genome,results_ncRNA_genome)
#   write.csv(results_noncoding_genome,paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricNoncoding/",trait,"_noncoding_sig.csv"),row.names = FALSE)
# }
# 
# 
# 
