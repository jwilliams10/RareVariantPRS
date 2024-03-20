rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_ncRNA_es_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/RareVariant_PRS/Merge_ncRNA_es_Binary.sh -icmd="bash Merge_ncRNA_es_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/GeneCentric_ncRNA/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2

Asthma_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  Asthma_ncRNA_coding_es <- rbind(Asthma_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/Asthma_ncRNA_sig_array",i,".csv")))
}

T2D_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  T2D_ncRNA_coding_es <- rbind(T2D_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/T2D_ncRNA_sig_array",i,".csv")))
}

Breast_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  Breast_ncRNA_coding_es <- rbind(Breast_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/Breast_ncRNA_sig_array",i,".csv")))
}

Prostate_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  Prostate_ncRNA_coding_es <- rbind(Prostate_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/Prostate_ncRNA_sig_array",i,".csv")))
}

CAD_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  CAD_ncRNA_coding_es <- rbind(CAD_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/CAD_ncRNA_sig_array",i,".csv")))
}

write.csv(Asthma_ncRNA_coding_es,row.names = FALSE,file = paste0("Asthma_ncRNA_sig_es.csv"))
write.csv(T2D_ncRNA_coding_es,row.names = FALSE,file = paste0("T2D_ncRNA_sig_es.csv"))
write.csv(Breast_ncRNA_coding_es,row.names = FALSE,file = paste0("Breast_ncRNA_sig_es.csv"))
write.csv(Prostate_ncRNA_coding_es,row.names = FALSE,file = paste0("Prostate_ncRNA_sig_es.csv"))
write.csv(CAD_ncRNA_coding_es,row.names = FALSE,file = paste0("CAD_ncRNA_sig_es.csv"))

system("rm -r GeneCentric_ncRNA_es/")
