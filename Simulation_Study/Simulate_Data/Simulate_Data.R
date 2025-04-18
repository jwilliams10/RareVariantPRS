rm(list = ls())
set.seed(1340)
gc()

library(bigsnpr)
library(dplyr)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

if(file.exists("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")){
  system("rm /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")
  system("rm /data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.bk")
}

rds <- bigsnpr::snp_readBed("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.bed") 

# Loading the data from backing files
common_variants <- snp_attach("/gpfs/gsfs12/users/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.rds")

fam_file <- read.table("/data/williamsjacr/UKB_WES_Simulation/chr22_fulldata/chr22_filtered_common.fam", quote="\"", comment.char="")

number_snps <- dim(common_variants$genotypes)[2]

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

causalprop_vec <- c(0.2,0.05,0.01,0.001,0.0005)
scale <- c(0,1)

Y <- list()

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

h2_overall_EUR <- vector()
h2_overall_EAS <- vector()
h2_overall_SAS <- vector()
h2_overall_AMR <- vector()
h2_overall_AFR <- vector()

h2_common_EUR <- vector()
h2_common_EAS <- vector()
h2_common_SAS <- vector()
h2_common_AMR <- vector()
h2_common_AFR <- vector()

MAF_common_EUR <- vector()
MAF_common_EAS <- vector()
MAF_common_SAS <- vector()
MAF_common_AMR <- vector()
MAF_common_AFR <- vector()

Var_common_EUR <- vector()
Var_common_EAS <- vector()
Var_common_SAS <- vector()
Var_common_AMR <- vector()
Var_common_AFR <- vector()

Average_Burden_rare_EUR <- vector()
Average_Burden_rare_EAS <- vector()
Average_Burden_rare_SAS <- vector()
Average_Burden_rare_AMR <- vector()
Average_Burden_rare_AFR <- vector()

Var_Burden_rare_EUR <- vector()
Var_Burden_rare_EAS <- vector()
Var_Burden_rare_SAS <- vector()
Var_Burden_rare_AMR <- vector()
Var_Burden_rare_AFR <- vector()

h2_rare_EUR <- vector()
h2_rare_EAS <- vector()
h2_rare_SAS <- vector()
h2_rare_AMR <- vector()
h2_rare_AFR <- vector()

count <- 1

h2_common <- 0.05
h2_rare <- 0.05/4

for(j in 1:length(causalprop_vec)){
  number_causal_snps <- round(causalprop_vec[j]*number_snps)
  causal_snps <- sample(1:number_snps,number_causal_snps)
  
  beta_snps <- rnorm(number_causal_snps,mean = 0,sqrt(h2_common/number_causal_snps))
  
  g_snps <- common_variants$genotypes[,causal_snps]
  for(i in 1:ncol(g_snps)){
    g_snps[is.na(g_snps[,i]),i] <- 0
  }
  scaled_causal_snps <- scale(g_snps)
  
  number_causal_genes <- round(causalprop_vec[j]*number_genes)
  causal_genes <- sample(1:number_genes,number_causal_genes)
  
  scaled_causal_genes <- scale(G_star_gene_centric_coding[,causal_genes,drop = FALSE])
  
  beta_genes <- rnorm(number_causal_genes,mean = 0,sqrt(h2_rare/number_causal_genes))
  
  for(q in 1:length(scale)){
    for(l in 1:100){
      
      
      if(scale[q] == 0){
        # Var(G Beta) = c c/v -> 0.05
        v <- var(g_snps%*%matrix(beta_snps,ncol = 1))/h2_common
        Y_hat_1 <- data.frame(IDs = fam_file[,2],Y_hat_Common = g_snps%*%matrix(beta_snps/sqrt(as.numeric(v)),ncol = 1))
        
      }else{
        v <- var(scaled_causal_snps%*%matrix(beta_snps,ncol = 1))/h2_common
        Y_hat_1 <- data.frame(IDs = fam_file[,2],Y_hat_Common = scaled_causal_snps%*%matrix(beta_snps/sqrt(as.numeric(v)),ncol = 1)) 
      }
      
      if(scale[q] == 0){
        # Var(G Beta) = c c/v -> 0.05
        v <- var(G_star_gene_centric_coding[,causal_genes,drop = FALSE]%*%matrix(beta_genes,ncol = 1))/h2_rare
        Y_hat_2 <- data.frame(IDs = as.numeric(phenotype.id.merge[,1]),Y_hat_Rare = G_star_gene_centric_coding[,causal_genes,drop = FALSE]%*%matrix(beta_genes/sqrt(as.numeric(v)),ncol = 1))
      }else{
        v <- var(scaled_causal_genes%*%matrix(beta_genes,ncol = 1))/h2_rare
        Y_hat_2 <- data.frame(IDs = as.numeric(phenotype.id.merge[,1]),Y_hat_Rare = scaled_causal_genes%*%matrix(beta_genes/sqrt(as.numeric(v)),ncol = 1)) 
      }
    
      ### Combine
      
      Y_hat <- inner_join(Y_hat_1,Y_hat_2)
      
      epsilon <- rnorm(nrow(Y_hat),mean = 0,sd = sqrt(1 - h2_common - h2_rare))
      
      Y_raw <- Y_hat$Y_hat_Common + Y_hat$Y_hat_Rare
      
      Y_tmp <- Y_raw + epsilon
      
      h2_common_EUR[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])
      h2_common_SAS[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])
      h2_common_AFR[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])
      h2_common_AMR[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])
      h2_common_EAS[count] <- var(Y_hat$Y_hat_Common[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])
      
      h2_rare_EUR[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])
      h2_rare_SAS[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])
      h2_rare_AFR[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])
      h2_rare_AMR[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])
      h2_rare_EAS[count] <- var(Y_hat$Y_hat_Rare[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])
      
      MAF_common_EUR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_SAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_AFR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_AMR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],],2,function(x){sum(x)/(2*length(x))}))
      MAF_common_EAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],],2,function(x){sum(x)/(2*length(x))}))
      
      Var_common_EUR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],],2,var))
      Var_common_SAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],],2,var))
      Var_common_AFR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],],2,var))
      Var_common_AMR[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],],2,var))
      Var_common_EAS[count] <- mean(apply(g_snps[Y_hat_1$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],],2,var))
      
      Average_Burden_rare_EUR[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],causal_genes,drop = FALSE]))
      Average_Burden_rare_SAS[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],causal_genes,drop = FALSE]))
      Average_Burden_rare_AFR[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],causal_genes,drop = FALSE]))
      Average_Burden_rare_AMR[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],causal_genes,drop = FALSE]))
      Average_Burden_rare_EAS[count] <- mean(colMeans(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],causal_genes,drop = FALSE]))
      
      Var_Burden_rare_EUR[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],causal_genes,drop = FALSE],2,var))
      Var_Burden_rare_SAS[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],causal_genes,drop = FALSE],2,var))
      Var_Burden_rare_AFR[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],causal_genes,drop = FALSE],2,var))
      Var_Burden_rare_AMR[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],causal_genes,drop = FALSE],2,var))
      Var_Burden_rare_EAS[count] <- mean(apply(G_star_gene_centric_coding[Y_hat_2$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],causal_genes,drop = FALSE],2,var))
      
      h2_overall_EUR[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"]])
      h2_overall_SAS[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"]])
      h2_overall_AFR[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"]])
      h2_overall_AMR[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"]])
      h2_overall_EAS[count] <- var(Y_raw[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])/var(Y_tmp[Y_hat$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"]])
      
      Y[[count]] <- data.frame(IDs = Y_hat$IDs,Y=Y_tmp)
      
      print(count)
      
      count <- count + 1
    }
  }
}

h2_dat <- data.frame(causal_prop = rep(rep(causalprop_vec,each = 200),5),scaled = rep(rep(rep(c("No","Yes"),each = 100),5),5),
  Ancestry = rep(c("EUR","SAS","AFR","AMR","EAS"),each = 1000),
  h2_common = c(h2_common_EUR,h2_common_SAS,h2_common_AFR,h2_common_AMR,h2_common_EAS),
  h2_rare = c(h2_rare_EUR,h2_rare_SAS,h2_rare_AFR,h2_rare_AMR,h2_rare_EAS),
  h2_overall = c(h2_overall_EUR,h2_overall_SAS,h2_overall_AFR,h2_overall_AMR,h2_overall_EAS),
  MAF_common = c(MAF_common_EUR,MAF_common_SAS,MAF_common_AFR,MAF_common_AMR,MAF_common_EAS),
  Var_common = c(Var_common_EUR,Var_common_SAS,Var_common_AFR,Var_common_AMR,Var_common_EAS),
  Average_Burden_rare = c(Average_Burden_rare_EUR,Average_Burden_rare_SAS,Average_Burden_rare_AFR,Average_Burden_rare_AMR,Average_Burden_rare_EAS),
  Var_Burden_rare = c(Var_Burden_rare_EUR,Var_Burden_rare_SAS,Var_Burden_rare_AFR,Var_Burden_rare_AMR,Var_Burden_rare_EAS))

write.csv(h2_dat,file = "/data/williamsjacr/UKB_WES_Simulation/Sim_Characteristics.csv",row.names = FALSE)

lapply(Y,function(x){var(x$Y)})
summary(unlist(lapply(Y,function(x){var(x$Y)})))

var(Y[[1]]$Y)
var(Y[[i]]$Y)

save(Y,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/Y_n_",nrow(Y_hat),"_h2_common_",h2_common,"_h2_rare_",h2_rare,".RData"))
save(Y,file = paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/Y_n_",nrow(Y_hat),"_h2_common_",h2_common,"_h2_rare_",h2_rare,".RData"))
