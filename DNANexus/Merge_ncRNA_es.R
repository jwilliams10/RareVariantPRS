rm(list = ls())

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Merge_ncRNA_es.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Continuous/RareVariant_PRS/Merge_ncRNA_es.sh -icmd="bash Merge_ncRNA_es.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/GeneCentric_ncRNA/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x2

BMI_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  BMI_ncRNA_coding_es <- rbind(BMI_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/BMI_ncRNA_sig_array",i,".csv")))
}

LDL_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  LDL_ncRNA_coding_es <- rbind(LDL_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/LDL_ncRNA_sig_array",i,".csv")))
}

HDL_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  HDL_ncRNA_coding_es <- rbind(HDL_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/HDL_ncRNA_sig_array",i,".csv")))
}

TC_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  TC_ncRNA_coding_es <- rbind(TC_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/TC_ncRNA_sig_array",i,".csv")))
}

logTG_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  logTG_ncRNA_coding_es <- rbind(logTG_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/logTG_ncRNA_sig_array",i,".csv")))
}

Height_ncRNA_coding_es <- NULL
for(i in c(1:200)){
  Height_ncRNA_coding_es <- rbind(Height_ncRNA_coding_es,read.csv(paste0("GeneCentric_ncRNA_es/Height_ncRNA_sig_array",i,".csv")))
}

write.csv(BMI_ncRNA_coding_es,row.names = FALSE,file = paste0("BMI_ncRNA_sig_es.csv"))
write.csv(LDL_ncRNA_coding_es,row.names = FALSE,file = paste0("LDL_ncRNA_sig_es.csv"))
write.csv(HDL_ncRNA_coding_es,row.names = FALSE,file = paste0("HDL_ncRNA_sig_es.csv"))
write.csv(TC_ncRNA_coding_es,row.names = FALSE,file = paste0("TC_ncRNA_sig_es.csv"))
write.csv(logTG_ncRNA_coding_es,row.names = FALSE,file = paste0("logTG_ncRNA_sig_es.csv"))
write.csv(Height_ncRNA_coding_es,row.names = FALSE,file = paste0("Height_ncRNA_sig_es.csv"))

system("rm -r GeneCentric_ncRNA_es/")
