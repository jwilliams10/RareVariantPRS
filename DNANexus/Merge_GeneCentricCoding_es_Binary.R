rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_GeneCentricCoding_es_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_GeneCentricCoding_es_Binary.sh -icmd="bash Merge_GeneCentricCoding_es_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentricCoding/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2


Asthma_coding_es <- NULL
for(i in c(1:400)){
  Asthma_coding_es <- rbind(Asthma_coding_es,read.csv(paste0("GeneCentricCoding_es/Asthma_coding_sig_array",i,".csv")))
}

T2D_coding_es <- NULL
for(i in c(1:400)){
  T2D_coding_es <- rbind(T2D_coding_es,read.csv(paste0("GeneCentricCoding_es/T2D_coding_sig_array",i,".csv")))
}

Breast_coding_es <- NULL
for(i in c(1:400)){
  Breast_coding_es <- rbind(Breast_coding_es,read.csv(paste0("GeneCentricCoding_es/Breast_coding_sig_array",i,".csv")))
}

Prostate_coding_es <- NULL
for(i in c(1:400)){
  Prostate_coding_es <- rbind(Prostate_coding_es,read.csv(paste0("GeneCentricCoding_es/Prostate_coding_sig_array",i,".csv")))
}

CAD_coding_es <- NULL
for(i in c(1:400)){
  CAD_coding_es <- rbind(CAD_coding_es,read.csv(paste0("GeneCentricCoding_es/CAD_coding_sig_array",i,".csv")))
}

write.csv(Asthma_coding_es,row.names = FALSE,file = paste0("Asthma_coding_sig_es.csv"))
write.csv(T2D_coding_es,row.names = FALSE,file = paste0("T2D_coding_sig_es.csv"))
write.csv(Breast_coding_es,row.names = FALSE,file = paste0("Breast_coding_sig_es.csv"))
write.csv(Prostate_coding_es,row.names = FALSE,file = paste0("Prostate_coding_sig_es.csv"))
write.csv(CAD_coding_es,row.names = FALSE,file = paste0("CAD_coding_sig_es.csv"))

system("rm -r GeneCentricCoding_es/")
