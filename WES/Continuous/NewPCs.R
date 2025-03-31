rm(list = ls())
library(dplyr)

All_Train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")

All_Train <- subset(All_Train,select = -c(pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10))
All_Train_New <- read.delim("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/All_Train.txt")
All_Tune_New <- read.delim("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/All_Tune.txt")
All_Validation_New <- read.delim("/data/williamsjacr/RICE_CVx_WES/RICE_CVx/All_Validation.txt")
All_New <- rbind(All_Train_New,All_Tune_New,All_Validation_New)

All_Train <- inner_join(All_Train,All_New[c("IID",paste0("pc",1:10))])

write.table(All_Train,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Train_NewPCs.txt",sep = '\t',row.names = FALSE,quote = FALSE)
