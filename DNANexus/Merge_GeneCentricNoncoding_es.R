rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Merge_GeneCentricNoncoding_es.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Merge_GeneCentricNoncoding_es.sh -icmd="bash Merge_GeneCentricNoncoding_es.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricNoncoding/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2

BMI_noncoding_es <- NULL
for(i in c(1:1600)){
  BMI_noncoding_es <- rbind(BMI_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/BMI_noncoding_sig_array",i,".csv")))
}

LDL_noncoding_es <- NULL
for(i in c(1:1600)){
  LDL_noncoding_es <- rbind(LDL_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/LDL_noncoding_sig_array",i,".csv")))
}

HDL_noncoding_es <- NULL
for(i in c(1:1600)){
  HDL_noncoding_es <- rbind(HDL_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/HDL_noncoding_sig_array",i,".csv")))
}

TC_noncoding_es <- NULL
for(i in c(1:1600)){
  TC_noncoding_es <- rbind(TC_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/TC_noncoding_sig_array",i,".csv")))
}

logTG_noncoding_es <- NULL
for(i in c(1:1600)){
  logTG_noncoding_es <- rbind(logTG_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/logTG_noncoding_sig_array",i,".csv")))
}

Height_noncoding_es <- NULL
for(i in c(1:1600)){
  Height_noncoding_es <- rbind(Height_noncoding_es,read.csv(paste0("GeneCentricNoncoding_es/Height_noncoding_sig_array",i,".csv")))
}

write.csv(BMI_noncoding_es,row.names = FALSE,file = paste0("BMI_noncoding_sig_sig_es.csv"))
write.csv(LDL_noncoding_es,row.names = FALSE,file = paste0("LDL_noncoding_sig_es.csv"))
write.csv(HDL_noncoding_es,row.names = FALSE,file = paste0("HDL_noncoding_sig_es.csv"))
write.csv(TC_noncoding_es,row.names = FALSE,file = paste0("TC_noncoding_sig_es.csv"))
write.csv(logTG_noncoding_es,row.names = FALSE,file = paste0("logTG_noncoding_sig_es.csv"))
write.csv(Height_noncoding_es,row.names = FALSE,file = paste0("Height_noncoding_sig_es.csv"))

system("rm -r GeneCentricNoncoding_es/")
