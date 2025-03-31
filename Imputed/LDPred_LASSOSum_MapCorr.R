rm(list=ls())
time <- system.time({
  library(data.table)
  library(dplyr)
  library(pROC)
  library(bigsnpr)
  library(bigsparser)
  library(readr)
  library(boot)
  library(bigstatsr)
  
  snp_subset <- function(x, ind.row = bigstatsr::rows_along(x$genotypes), ind.col = bigstatsr::cols_along(x$genotypes), backingfile = NULL){
    G <- x$genotypes
    ind.row <- bigstatsr::rows_along(G)[ind.row]
    ind.col <- bigstatsr::cols_along(G)[ind.col]
    # check_args()
    if (is.null(x$fam)) {
      new_fam <- NULL
    }
    else {
      new_fam <- x$fam[ind.row, , drop = FALSE]
      rownames(new_fam) <- bigstatsr::rows_along(new_fam)
    }
    if (is.null(x$map)) {
      new_map <- NULL
    }
    else {
      new_map <- x$map[ind.col, , drop = FALSE]
      rownames(new_map) <- bigstatsr::rows_along(new_map)
    }
    if (is.null(backingfile)) 
      backingfile <- bigstatsr:::getNewFile(x, "sub")
    G2 <- bigstatsr:::FBM.code256(nrow = length(ind.row), ncol = length(ind.col), 
                                  code = G$code256, init = NULL, backingfile = backingfile, 
                                  create_bk = TRUE)
    bigsnpr:::replaceSNP(G2, G, rowInd = ind.row, colInd = ind.col)
    snp.list <- structure(list(genotypes = G2, fam = new_fam, 
                               map = new_map), class = "bigSNP")
    rds <- bigstatsr:::sub_bk(G2$backingfile, ".rds")
    saveRDS(snp.list, rds)
    rds
  }
  
  fill <- "all_chr_EUR_reference"
  
  ldr <- 3/1000
  ncores <- 1
  
  map <- NULL
  
  obj.bigSNP <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/BEDFiles/",fill,".rds"))
  NCORES <-  1
  
  EUR_LDBlocks <- read.csv("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/EUR_LDBlocks.txt", sep="")
  subset_logic_list <- lapply(1:nrow(EUR_LDBlocks),function(x){return(which((obj.bigSNP$map$physical.pos >= EUR_LDBlocks$start[x]) & (obj.bigSNP$map$physical.pos <= EUR_LDBlocks$stop[x]) & (obj.bigSNP$map$chromosome == EUR_LDBlocks$chr[x])))})
  empty_lists <- which(unlist(lapply(subset_logic_list, length)) == 0)
  subset_logic_list <- subset_logic_list[-empty_lists]
  
  for(i in 1:length(subset_logic_list)){
    snp_subset(obj.bigSNP,ind.row = 1:3000,ind.col = subset_logic_list[[i]],backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i)) 
    
    obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".rds"))
    
    map_new <- obj.bigSNP_new$map[-c(3)]
    names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
    
    G   <- obj.bigSNP_new$genotypes
    CHR <- obj.bigSNP_new$map$chromosome
    POS <- obj.bigSNP_new$map$physical.pos
    POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDPred2_Genetic_Mappings/", ncores = ncores)
    
    corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
    
    if(anyNA(corr0@x)){
      b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
      b <- as.numeric(names(table(b))[table(b) > 1])
      
      file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".rds"))
      file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".bk"))
      
      snp_subset(obj.bigSNP_new,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i))
      
      obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".rds"))
      
      map_new <- obj.bigSNP_new$map[-c(3)]
      names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
      
      G   <- obj.bigSNP_new$genotypes
      CHR <- obj.bigSNP_new$map$chromosome
      POS <- obj.bigSNP_new$map$physical.pos
      POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDPred2_Genetic_Mappings/", ncores = ncores)
      
      corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
      while(anyNA(corr0@x)){
        b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
        b <- as.numeric(names(table(b))[table(b) > 2])
        
        if(length(b) == 0){
          b <- Matrix::which(is.nan(corr0), arr.ind = TRUE)
          b <- as.numeric(names(table(b))[table(b) > 1])
        }
        
        file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".rds"))
        file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".bk"))
        
        snp_subset(obj.bigSNP_new,ind.row = 1:3000,ind.col = -b,backingfile = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i))
        
        obj.bigSNP_new <- snp_attach(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".rds"))
        
        map_new <- obj.bigSNP_new$map[-c(3)]
        names(map_new) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")
        
        G   <- obj.bigSNP_new$genotypes
        CHR <- obj.bigSNP_new$map$chromosome
        POS <- obj.bigSNP_new$map$physical.pos
        POS2 <- snp_asGeneticPos(CHR, POS, dir ="/data/williamsjacr/UKB_WES_Phenotypes/Continuous/LDPred2_Genetic_Mappings/", ncores = ncores)
        
        corr0 <- snp_cor(G,infos.pos = POS2, size =  ldr)
      }
    }
    
    file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".rds"))
    file.remove(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/LDpred2_LDFiles/tmp_",fill,"_",i,".bk"))
    
    map <- rbind(map, map_new)
    
    if(i == 1){
      corr <- as_SFBM(corr0, paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/",fill,"_Corr_Backingfile"), compact = TRUE)
    }else{
      corr$add_columns(corr0, nrow(corr))
    }
  }
  save(map,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/",fill,"_Map.RData"))
  save(corr,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/",fill,"_Corr.RData"))
})[3]
save(time,file = "/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2_Corr_Time.RData")