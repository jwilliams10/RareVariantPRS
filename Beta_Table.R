rm(list = ls())
WES_Results_Binary <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Binary.csv")
WES_Results_Continuous <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Continuous.csv")

full_results <- rbind(WES_Results_Binary,WES_Results_Continuous)

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","LASSOSum","LDPred","CV","RV"),]
full_results$Method[full_results$Method == "CV"] <- "RICE-CV" 
full_results$Method[full_results$Method == "RV"] <- "RICE-RV" 
full_results$Method[full_results$Method == "LDPred"] <- "LDpred2"
full_results$Method[full_results$Method == "LASSOSum"] <- "Lassosum2"

full_results$trait <- factor(full_results$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D","BMI","Height","HDL","LDL","logTG","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-CV","RICE-RV"))

full_results <- full_results[order(full_results$trait,full_results$ancestry,full_results$Method),]

write.csv(full_results,file = "Downloads/WES_Results.csv",row.names = FALSE)

WGS_Results_Binary <- read.csv("~/Desktop/RareVariantPRS_Results/WGS_Results_Binary.csv")
WGS_Results_Continuous <- read.csv("~/Desktop/RareVariantPRS_Results/WGS_Results_Continuous.csv")

full_results <- rbind(WGS_Results_Binary,WGS_Results_Continuous)

full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","LASSOSum","LDPred","CV","RV"),]
full_results$Method[full_results$Method == "CV"] <- "RICE-CV" 
full_results$Method[full_results$Method == "RV"] <- "RICE-RV" 
full_results$Method[full_results$Method == "LDPred"] <- "LDpred2"
full_results$Method[full_results$Method == "LASSOSum"] <- "Lassosum2"

full_results$trait <- factor(full_results$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D","BMI","Height","HDL","LDL","logTG","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-CV","RICE-RV"))

full_results <- full_results[order(full_results$trait,full_results$ancestry,full_results$Method),]

write.csv(full_results,file = "Downloads/WGS_Results.csv",row.names = FALSE)

AoU_Results <- read.csv("~/Desktop/RareVariantPRS_Results/AoU_Results.csv")

full_results <- AoU_Results

full_results <- full_results[full_results$Method %in% c("CTSLEB","JointPRS","PROSPER","RICE-CV","RICE-RV"),]
full_results$Method[full_results$Method == "CTSLEB"] <- "CT-SLEB"

full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","logTG","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EAS","EUR","MID","SAS"))
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-CV","RICE-RV"))

full_results <- full_results[order(full_results$trait,full_results$ancestry,full_results$Method),]

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

write.csv(full_results,file = "Downloads/AoU_Results.csv",row.names = FALSE)