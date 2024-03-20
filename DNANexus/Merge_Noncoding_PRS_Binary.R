rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_Noncoding_PRS_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_Noncoding_PRS_Binary.sh -icmd="bash Merge_Noncoding_PRS_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricNoncoding/ --priority high --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x8

tmp_info <- NULL

for(i in 1:5884){
  arrayid <- i
  arrayid_original <- arrayid
  
  if(arrayid>2942){
    Burden <- 1
    arrayid <- arrayid - 2942
  }else{
    Burden <- 0
  }
  
  if(arrayid <= 22){
    threshold <- 1
  }else if(arrayid <= 44){
    threshold <- 2
  }else if(arrayid <= 66){
    threshold <- 3
  }else if(arrayid <= 88){
    threshold <- 4
  }else if(arrayid <= 110){
    threshold <- 5
  }else if(arrayid <= 132){
    threshold <- 6
  }else if(arrayid <= 154){
    threshold <- 7
  }else if(arrayid <= 176){
    threshold <- 8
  }else if(arrayid <= 198){
    threshold <- 9
  }else if(arrayid <= 220){
    threshold <- 10
  }else if(arrayid <= 242){
    threshold <- 11
  }else if(arrayid <= 342){
    threshold <- 12
  }else if(arrayid <= 542){
    threshold <- 13
  }else if(arrayid <= 1342){
    threshold <- 14
  }else{
    threshold <- 15
  }
  tmp_info <- rbind(tmp_info,data.frame(arrayid = arrayid_original,Burden = Burden,Threshold = threshold))
}

rm(list=setdiff(ls(), "tmp_info"))

for(trait in c("Asthma","T2D","Breast","Prostate","CAD")){
  PRS <- read.csv(paste0("GeneCentricNoncoding_PRS/",trait,"_Noncoding_PRS5884.csv")) 
  IDs <- PRS$ID
  PRS_Final_Burden <- matrix(0,nrow = nrow(PRS),ncol = 15)
  PRS_Final_STAARO <- matrix(0,nrow = nrow(PRS),ncol = 15)
  
  for(array in 1:5884){
    Burden <- tmp_info[array,"Burden"]
    Threshold <- tmp_info[array,"Threshold"]
    
    tmp <- read.csv(paste0("GeneCentricNoncoding_PRS/",trait,"_Noncoding_PRS",array,".csv"))
    
    if(Burden == 1){
      PRS_Final_Burden[,Threshold] <- PRS_Final_Burden[,Threshold] + tmp$PRS
    }else{
      PRS_Final_STAARO[,Threshold] <- PRS_Final_STAARO[,Threshold] + tmp$PRS
    }
  }
  PRS_Final_Burden <- data.frame(IID = IDs, PRS = PRS_Final_Burden)
  colnames(PRS_Final_Burden) <- c("IID",paste0("Burden_Threshold",1:15))
  
  PRS_Final_STAARO <- data.frame(IID = IDs, PRS = PRS_Final_STAARO)
  colnames(PRS_Final_STAARO) <- c("IID",paste0("STAARO_Threshold",1:15))
  
  write.csv(PRS_Final_Burden,file = paste0(trait,"_Noncoding_Burden_PRS.csv"),row.names = FALSE)
  write.csv(PRS_Final_STAARO,file = paste0(trait,"_Noncoding_STAARO_PRS.csv"),row.names = FALSE)
}

system("rm -r GeneCentricNoncoding_PRS/")
