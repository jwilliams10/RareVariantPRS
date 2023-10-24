rm(list = ls())
set.seed(1330)
gc()

library(bigsnpr)
library(dplyr)

#### Common

if(file.exists("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")){
  system("rm /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")
  system("rm /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.bk")
}

rds <- bigsnpr::snp_readBed("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.bed") 

# Loading the data from backing files
common_variants <- snp_attach("/gpfs/gsfs12/users/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")

fam_file <- read.delim("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam", header=FALSE)

number_snps <- dim(common_variants$genotypes)[2]
number_causal_snps <- round(0.01*number_snps)

h2_common <- 0.05

causal_snps <- sample(1:number_snps,number_causal_snps)

g_snps <- common_variants$genotypes[,causal_snps]
for(i in 1:ncol(g_snps)){
  g_snps[is.na(g_snps[,i]),i] <- 0
}
scaled_causal_snps <- scale(g_snps)

beta_snps <- rnorm(number_causal_snps,mean = 0,sqrt(h2_common/number_causal_snps))

Y_hat_1 <- data.frame(IDs = fam_file[,1],Y_hat_Common = scaled_causal_snps%*%matrix(beta_snps,ncol = 1))

### Rare

load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/g_star/GeneCentric_Coding.RData")

obj_nullmodel <- get(load("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/samplenull_model.RData"))

genofile <- seqOpen("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/gds/full_gds22.gds")

id.genotype <- seqGetData(genofile,"sample.id")
id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
phenotype.id <- as.character(obj_nullmodel$id_include)
phenotype.id.merge <- data.frame(phenotype.id)
phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))

seqClose(genofile)

number_genes <- dim(G_star_gene_centric_coding)[2]
number_causal_genes <- round(0.01*number_genes)

h2_rare <- 0.05/4

causal_genes <- sample(1:number_genes,number_causal_genes)

scaled_causal_genes <- scale(G_star_gene_centric_coding[,causal_genes])

beta_genes <- rnorm(number_causal_genes,mean = 0,sqrt(h2_rare/number_causal_genes))

Y_hat_2 <- data.frame(IDs = as.numeric(phenotype.id.merge[,1]),Y_hat_Rare = scaled_causal_genes%*%matrix(beta_genes,ncol = 1))

### Combine

Y_hat <- inner_join(Y_hat_1,Y_hat_2)

Y <- list()

for(i in 1:10){
  epsilon <- rnorm(nrow(Y_hat),mean = 0,sd = sqrt(1 - h2_common - h2_rare))
  
  Y_tmp <- Y_hat$Y_hat_Common + Y_hat$Y_hat_Rare + epsilon
  
  Y[[i]] <- data.frame(IDs = Y_hat$IDs,Y=Y_tmp)
}

var(Y[[1]]$Y)
var(Y[[i]]$Y)

save(Y,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/Y_n_",nrow(Y_hat),"_h2_common_",h2_common,"_h2_rare_",h2_rare,".RData"))
