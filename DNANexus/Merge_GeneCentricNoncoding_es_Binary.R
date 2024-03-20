rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_GeneCentricNoncoding_es_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_GeneCentricNoncoding_es_Binary.sh -icmd="bash Merge_GeneCentricNoncoding_es_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricNoncoding/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2

Asthma_noncoding_es <- NULL
for(i in c(1:1600)){
  Asthma_noncoding_es <- rbind(Asthma_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/Asthma_noncoding_sig_array",i,".csv")))
}

T2D_noncoding_es <- NULL
for(i in c(1:1600)){
  T2D_noncoding_es <- rbind(T2D_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/T2D_noncoding_sig_array",i,".csv")))
}

Breast_noncoding_es <- NULL
for(i in c(1:1600)){
  Breast_noncoding_es <- rbind(Breast_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/Breast_noncoding_sig_array",i,".csv")))
}

Prostate_noncoding_es <- NULL
for(i in c(1:1600)){
  Prostate_noncoding_es <- rbind(Prostate_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/Prostate_noncoding_sig_array",i,".csv")))
}

CAD_noncoding_es <- NULL
for(i in c(1:1600)){
  CAD_noncoding_es <- rbind(CAD_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/CAD_noncoding_sig_array",i,".csv")))
}

write.csv(Asthma_noncoding_es,row.names = FALSE,file = paste0("Asthma_noncoding_sig_es.csv"))
write.csv(T2D_noncoding_es,row.names = FALSE,file = paste0("T2D_noncoding_sig_es.csv"))
write.csv(Breast_noncoding_es,row.names = FALSE,file = paste0("Breast_noncoding_sig_es.csv"))
write.csv(Prostate_noncoding_es,row.names = FALSE,file = paste0("Prostate_noncoding_sig_es.csv"))
write.csv(CAD_noncoding_es,row.names = FALSE,file = paste0("CAD_noncoding_sig_es.csv"))

system("rm -r GeneCentricNoncoding_es/")
