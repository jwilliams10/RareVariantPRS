rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Merge_GeneCentricCoding_es.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Merge_GeneCentricCoding_es.sh -icmd="bash Merge_GeneCentricCoding_es.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentricCoding/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2

BMI_coding_es <- NULL
for(i in c(1:400)){
  BMI_coding_es <- rbind(BMI_coding_es,read.csv(paste0("GeneCentricCoding_es/BMI_coding_sig_array",i,".csv")))
}

LDL_coding_es <- NULL
for(i in c(1:400)){
  LDL_coding_es <- rbind(LDL_coding_es,read.csv(paste0("GeneCentricCoding_es/LDL_coding_sig_array",i,".csv")))
}

HDL_coding_es <- NULL
for(i in c(1:400)){
  HDL_coding_es <- rbind(HDL_coding_es,read.csv(paste0("GeneCentricCoding_es/HDL_coding_sig_array",i,".csv")))
}

TC_coding_es <- NULL
for(i in c(1:400)){
  TC_coding_es <- rbind(TC_coding_es,read.csv(paste0("GeneCentricCoding_es/TC_coding_sig_array",i,".csv")))
}

logTG_coding_es <- NULL
for(i in c(1:400)){
  logTG_coding_es <- rbind(logTG_coding_es,read.csv(paste0("GeneCentricCoding_es/logTG_coding_sig_array",i,".csv")))
}

Height_coding_es <- NULL
for(i in c(1:400)){
  Height_coding_es <- rbind(Height_coding_es,read.csv(paste0("GeneCentricCoding_es/Height_coding_sig_array",i,".csv")))
}

write.csv(BMI_coding_es,row.names = FALSE,file = paste0("BMI_coding_sig_es.csv"))
write.csv(LDL_coding_es,row.names = FALSE,file = paste0("LDL_coding_sig_es.csv"))
write.csv(HDL_coding_es,row.names = FALSE,file = paste0("HDL_coding_sig_es.csv"))
write.csv(TC_coding_es,row.names = FALSE,file = paste0("TC_coding_sig_es.csv"))
write.csv(logTG_coding_es,row.names = FALSE,file = paste0("logTG_coding_sig_es.csv"))
write.csv(Height_coding_es,row.names = FALSE,file = paste0("Height_coding_sig_es.csv"))

system("rm -r GeneCentricCoding_es/")
