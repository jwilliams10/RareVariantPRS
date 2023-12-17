rm(list = ls())

Burden <- vector()
chr <- vector()
threshold <- vector()
for(i in 1:990){
  
  arrayid <- i
  
  if(arrayid>495){
    Burden[i] <- 1
    arrayid <- arrayid - 495
  }else{
    Burden[i] <- 0
  }
  
  if(arrayid <= 33){
    threshold[i] <- 1
  }else if(arrayid <= 66){
    threshold[i] <- 2
  }else if(arrayid <= 99){
    threshold[i] <- 3
  }else if(arrayid <= 132){
    threshold[i] <- 4
  }else if(arrayid <= 165){
    threshold[i] <- 5
  }else if(arrayid <= 198){
    threshold[i] <- 6
  }else if(arrayid <= 231){
    threshold[i] <- 7
  }else if(arrayid <= 264){
    threshold[i] <- 8
  }else if(arrayid <= 297){
    threshold[i] <- 9
  }else if(arrayid <= 330){
    threshold[i] <- 10
  }else if(arrayid <= 363){
    threshold[i] <- 11
  }else if(arrayid <= 396){
    threshold[i] <- 12
  }else if(arrayid <= 429){
    threshold[i] <- 13
  }else if(arrayid <= 462){
    threshold[i] <- 14
  }else{
    threshold[i] <- 15
  }
  
  chunk <- arrayid %% 33
  if(chunk == 0){chunk <- 33}
}

data_results <- data.frame(arrayid = 1:990,Burden = Burden, Threshold = threshold)

rm(list = setdiff(ls(),"data_results"))

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  obj_nullmodel_tune <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Tune_Null_Model.RData")))
  obj_nullmodel_validation <- get(load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/nullmodels_staar/",trait,"_Validation_Null_Model.RData")))
  
  full <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_330.csv"))
  
  idx_tune <- which(full$ID %in% obj_nullmodel_tune$id_include)
  idx_validation <- which(full$ID %in% obj_nullmodel_validation$id_include)
  
  burden <- 0
  threshold <- 1
  
  for(burden in c(0,1)){
    for(threshold in 1:15){
      
      arrayid <- data_results$arrayid[data_results$Burden == burden & data_results$Threshold == threshold]
      
      if(threshold == 1){
        if(burden == 0){
          for(a in arrayid){
            if(a == arrayid[1]){
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))
              STAARO_GeneCentric_Coding_Tune_PRS <- tmp[idx_tune,]
              STAARO_GeneCentric_Coding_Validation_PRS <- tmp[idx_validation,]
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))
              STAARO_GeneCentric_Noncoding_Tune_PRS <- tmp[idx_tune,]
              STAARO_GeneCentric_Noncoding_Validation_PRS <- tmp[idx_validation,]
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))
              STAARO_SlidingWindow_Tune_PRS <- tmp[idx_tune,]
              STAARO_SlidingWindow_Validation_PRS <- tmp[idx_validation,]
            }else{
              STAARO_GeneCentric_Coding_Tune_PRS[,2] <- STAARO_GeneCentric_Coding_Tune_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              STAARO_GeneCentric_Coding_Validation_PRS[,2] <- STAARO_GeneCentric_Coding_Validation_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              STAARO_GeneCentric_Noncoding_Tune_PRS[,2] <- STAARO_GeneCentric_Noncoding_Tune_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              STAARO_GeneCentric_Noncoding_Validation_PRS[,2] <- STAARO_GeneCentric_Noncoding_Validation_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              STAARO_SlidingWindow_Tune_PRS[,2] <- STAARO_SlidingWindow_Tune_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              STAARO_SlidingWindow_Validation_PRS[,2] <- STAARO_SlidingWindow_Validation_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
            }
          }
          
          STAARO_GeneCentric_Coding_Tune_PRS$ID <- full$ID[idx_tune]
          STAARO_GeneCentric_Coding_Validation_PRS$ID <- full$ID[idx_validation]
          
          colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c("IID","PRS_Threshold_1")
          colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c("IID","PRS_Threshold_1")
          
          STAARO_GeneCentric_Noncoding_Tune_PRS$ID <- full$ID[idx_tune]
          STAARO_GeneCentric_Noncoding_Validation_PRS$ID <- full$ID[idx_validation]
          
          colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c("IID","PRS_Threshold_1")
          colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c("IID","PRS_Threshold_1")
          
          STAARO_SlidingWindow_Tune_PRS$ID <- full$ID[idx_tune]
          STAARO_SlidingWindow_Validation_PRS$ID <- full$ID[idx_validation]
          
          colnames(STAARO_SlidingWindow_Tune_PRS) <- c("IID","PRS_Threshold_1")
          colnames(STAARO_SlidingWindow_Validation_PRS) <- c("IID","PRS_Threshold_1")
        }else{
          for(a in arrayid){
            if(a == arrayid[1]){
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))
              Burden_GeneCentric_Coding_Tune_PRS <- tmp[idx_tune,]
              Burden_GeneCentric_Coding_Validation_PRS <- tmp[idx_validation,]
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))
              Burden_GeneCentric_Noncoding_Tune_PRS <- tmp[idx_tune,]
              Burden_GeneCentric_Noncoding_Validation_PRS <- tmp[idx_validation,]
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))
              Burden_SlidingWindow_Tune_PRS <- tmp[idx_tune,]
              Burden_SlidingWindow_Validation_PRS <- tmp[idx_validation,]
            }else{
              Burden_GeneCentric_Coding_Tune_PRS[,2] <- Burden_GeneCentric_Coding_Tune_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              Burden_GeneCentric_Coding_Validation_PRS[,2] <- Burden_GeneCentric_Coding_Validation_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              Burden_GeneCentric_Noncoding_Tune_PRS[,2] <- Burden_GeneCentric_Noncoding_Tune_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              Burden_GeneCentric_Noncoding_Validation_PRS[,2] <- Burden_GeneCentric_Noncoding_Validation_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              Burden_SlidingWindow_Tune_PRS[,2] <- Burden_SlidingWindow_Tune_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              Burden_SlidingWindow_Validation_PRS[,2] <- Burden_SlidingWindow_Validation_PRS[,2] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
            }
          }
          Burden_GeneCentric_Coding_Tune_PRS$ID <- full$ID[idx_tune]
          Burden_GeneCentric_Coding_Validation_PRS$ID <- full$ID[idx_validation]
          
          colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c("IID","PRS_Threshold_1")
          colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c("IID","PRS_Threshold_1")
          
          Burden_GeneCentric_Noncoding_Tune_PRS$ID <- full$ID[idx_tune]
          Burden_GeneCentric_Noncoding_Validation_PRS$ID <- full$ID[idx_validation]
          
          colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c("IID","PRS_Threshold_1")
          colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c("IID","PRS_Threshold_1")
          
          Burden_SlidingWindow_Tune_PRS$ID <- full$ID[idx_tune]
          Burden_SlidingWindow_Validation_PRS$ID <- full$ID[idx_validation]
          
          colnames(Burden_SlidingWindow_Tune_PRS) <- c("IID","PRS_Threshold_1")
          colnames(Burden_SlidingWindow_Validation_PRS) <- c("IID","PRS_Threshold_1")
        }
      }else{
        if(burden == 0){
          for(a in arrayid){
            if(a == arrayid[1]){
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))
              STAARO_GeneCentric_Coding_Tune_PRS <- cbind(STAARO_GeneCentric_Coding_Tune_PRS,tmp[idx_tune,2])
              STAARO_GeneCentric_Coding_Validation_PRS <- cbind(STAARO_GeneCentric_Coding_Validation_PRS,tmp[idx_validation,2])
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))
              STAARO_GeneCentric_Noncoding_Tune_PRS <- cbind(STAARO_GeneCentric_Noncoding_Tune_PRS,tmp[idx_tune,2])
              STAARO_GeneCentric_Noncoding_Validation_PRS <- cbind(STAARO_GeneCentric_Noncoding_Validation_PRS,tmp[idx_validation,2])
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))
              STAARO_SlidingWindow_Tune_PRS <- cbind(STAARO_SlidingWindow_Tune_PRS,tmp[idx_tune,2])
              STAARO_SlidingWindow_Validation_PRS <- cbind(STAARO_SlidingWindow_Validation_PRS,tmp[idx_validation,2])
            }else{
              STAARO_GeneCentric_Coding_Tune_PRS[,threshold + 1] <- STAARO_GeneCentric_Coding_Tune_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              STAARO_GeneCentric_Coding_Validation_PRS[,threshold + 1] <- STAARO_GeneCentric_Coding_Validation_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              STAARO_GeneCentric_Noncoding_Tune_PRS[,threshold + 1] <- STAARO_GeneCentric_Noncoding_Tune_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              STAARO_GeneCentric_Noncoding_Validation_PRS[,threshold + 1] <- STAARO_GeneCentric_Noncoding_Validation_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              STAARO_SlidingWindow_Tune_PRS[,threshold + 1] <- STAARO_SlidingWindow_Tune_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              STAARO_SlidingWindow_Validation_PRS[,threshold + 1] <- STAARO_SlidingWindow_Validation_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
            }
          }
          colnames(STAARO_GeneCentric_Coding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Tune_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          colnames(STAARO_GeneCentric_Coding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Coding_Validation_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          
          colnames(STAARO_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Tune_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          colnames(STAARO_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(STAARO_GeneCentric_Noncoding_Validation_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          
          colnames(STAARO_SlidingWindow_Tune_PRS) <- c(colnames(STAARO_SlidingWindow_Tune_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          colnames(STAARO_SlidingWindow_Validation_PRS) <- c(colnames(STAARO_SlidingWindow_Validation_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
        }else{
          for(a in arrayid){
            if(a == arrayid[1]){
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))
              Burden_GeneCentric_Coding_Tune_PRS <- cbind(Burden_GeneCentric_Coding_Tune_PRS,tmp[idx_tune,2])
              Burden_GeneCentric_Coding_Validation_PRS <- cbind(Burden_GeneCentric_Coding_Validation_PRS,tmp[idx_validation,2])
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))
              Burden_GeneCentric_Noncoding_Tune_PRS <- cbind(Burden_GeneCentric_Noncoding_Tune_PRS,tmp[idx_tune,2])
              Burden_GeneCentric_Noncoding_Validation_PRS <- cbind(Burden_GeneCentric_Noncoding_Validation_PRS,tmp[idx_validation,2])
              
              tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))
              Burden_SlidingWindow_Tune_PRS <- cbind(Burden_SlidingWindow_Tune_PRS,tmp[idx_tune,2])
              Burden_SlidingWindow_Validation_PRS <- cbind(Burden_SlidingWindow_Validation_PRS,tmp[idx_validation,2])
            }else{
              Burden_GeneCentric_Coding_Tune_PRS[,threshold + 1] <- Burden_GeneCentric_Coding_Tune_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              Burden_GeneCentric_Coding_Validation_PRS[,threshold + 1] <- Burden_GeneCentric_Coding_Validation_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              Burden_GeneCentric_Noncoding_Tune_PRS[,threshold + 1] <- Burden_GeneCentric_Noncoding_Tune_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              Burden_GeneCentric_Noncoding_Validation_PRS[,threshold + 1] <- Burden_GeneCentric_Noncoding_Validation_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
              
              Burden_SlidingWindow_Tune_PRS[,threshold + 1] <- Burden_SlidingWindow_Tune_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_tune,2]
              Burden_SlidingWindow_Validation_PRS[,threshold + 1] <- Burden_SlidingWindow_Validation_PRS[,threshold + 1] + read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_PRS_Array_",a,".csv"))[idx_validation,2]
            }
          } 
          colnames(Burden_GeneCentric_Coding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Coding_Tune_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          colnames(Burden_GeneCentric_Coding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Coding_Validation_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          
          colnames(Burden_GeneCentric_Noncoding_Tune_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Tune_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          colnames(Burden_GeneCentric_Noncoding_Validation_PRS) <- c(colnames(Burden_GeneCentric_Noncoding_Validation_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          
          colnames(Burden_SlidingWindow_Tune_PRS) <- c(colnames(Burden_SlidingWindow_Tune_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
          colnames(Burden_SlidingWindow_Validation_PRS) <- c(colnames(Burden_SlidingWindow_Validation_PRS)[1:threshold],paste0("PRS_Threshold_",threshold))
        }
      }
    }
  }
  
  write.csv(STAARO_GeneCentric_Coding_Tune_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_STAARO_GeneCentric_Coding_Tune_PRS.csv"),row.names = FALSE)
  write.csv(STAARO_GeneCentric_Coding_Validation_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_STAARO_GeneCentric_Coding_Validation_PRS.csv"),row.names = FALSE)
  
  write.csv(STAARO_GeneCentric_Noncoding_Tune_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_STAARO_GeneCentric_Noncoding_Tune_PRS.csv"),row.names = FALSE)
  write.csv(STAARO_GeneCentric_Noncoding_Validation_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_STAARO_GeneCentric_Noncoding_Validation_PRS.csv"),row.names = FALSE)
  
  write.csv(STAARO_SlidingWindow_Tune_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_STAARO_SlidingWindow_Tune_PRS.csv"),row.names = FALSE)
  write.csv(STAARO_SlidingWindow_Validation_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_STAARO_SlidingWindow_Validation_PRS.csv"),row.names = FALSE)
  
  write.csv(Burden_GeneCentric_Coding_Tune_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_Burden_GeneCentric_Coding_Tune_PRS.csv"),row.names = FALSE)
  write.csv(Burden_GeneCentric_Coding_Validation_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_Burden_GeneCentric_Coding_Validation_PRS.csv"),row.names = FALSE)
  
  write.csv(Burden_GeneCentric_Noncoding_Tune_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_Burden_GeneCentric_Noncoding_Tune_PRS.csv"),row.names = FALSE)
  write.csv(Burden_GeneCentric_Noncoding_Validation_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricNoncoding/",trait,"_Burden_GeneCentric_Noncoding_Validation_PRS.csv"),row.names = FALSE)
  
  write.csv(Burden_SlidingWindow_Tune_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_Burden_SlidingWindow_Tune_PRS.csv"),row.names = FALSE)
  write.csv(Burden_SlidingWindow_Validation_PRS,file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/SlidingWindow/",trait,"_Burden_SlidingWindow_Validation_PRS.csv"),row.names = FALSE)
  
}
