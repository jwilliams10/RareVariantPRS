rm(list = ls())

All_Train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")

BMIadj_resid <- resid(lm(BMI~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = All_Train))
All_Train$BMIadj_norm[!is.na(All_Train$BMI)] <- qnorm((rank(BMIadj_resid,na.last="keep")-0.5)/length(BMIadj_resid))

Heightadj_resid <- resid(lm(Height~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = All_Train))
All_Train$Heightadj_norm[!is.na(All_Train$Height)] <- qnorm((rank(Heightadj_resid,na.last="keep")-0.5)/length(Heightadj_resid))

LDLadj_resid <- resid(lm(LDL~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = All_Train))
All_Train$LDLadj_norm[!is.na(All_Train$LDL)] <- qnorm((rank(LDLadj_resid,na.last="keep")-0.5)/length(LDLadj_resid))

HDLadj_resid <- resid(lm(HDL~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = All_Train))
All_Train$HDLadj_norm[!is.na(All_Train$HDL)] <- qnorm((rank(HDLadj_resid,na.last="keep")-0.5)/length(HDLadj_resid))

logTGadj_resid <- resid(lm(logTG~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = All_Train))
All_Train$logTGadj_norm[!is.na(All_Train$logTG)] <- qnorm((rank(logTGadj_resid,na.last="keep")-0.5)/length(logTGadj_resid))

TCadj_resid <- resid(lm(TC~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = All_Train))
All_Train$TCadj_norm[!is.na(All_Train$TC)] <- qnorm((rank(TCadj_resid,na.last="keep")-0.5)/length(TCadj_resid))

write.table(All_Train,file = "/data/williamsjacr/UKB_WES_Phenotypes/All_Train_RankNormal.txt",sep = '\t',row.names = FALSE,quote = FALSE)