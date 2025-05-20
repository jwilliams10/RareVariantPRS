rm(list = ls())

coding_qqplot_STAARB <- function(coding_sig,coding_sig_plof,coding_sig_plofds,coding_sig_missense,coding_sig_synonymous,coding_sig_disruptive_missense,trait){
  coding_sig$STAARB[coding_sig$STAARB == 0] <- 10^(-250)
  coding_minp <- min(coding_sig$STAARB)
  min_y <- ceiling(-log10(coding_minp)) + 1
  cex_point <- 1
  
  
  observed <- sort(coding_sig_plof$STAARB)
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  observed <- sort(coding_sig_plofds$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_missense$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### disruptive_missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_disruptive_missense$STAARB)
  
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(coding_sig_synonymous$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),ncol=1,bty="o",box.lwd=1,pch=0:4,cex=1.3,text.font=2)
  
}

print("WES Coding STAARB")

png(paste0("UKB_WES_Continuous_RV_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), 2, 6, byrow = TRUE))

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]
  
  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}

dev.off()

png(paste0("UKB_WES_Binary_RV_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2, 6, byrow = TRUE))

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]
  
  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}
dev.off()

rm(list = ls())

coding_qqplot_STAARB <- function(coding_sig,coding_sig_plof,coding_sig_plofds,coding_sig_missense,coding_sig_synonymous,coding_sig_disruptive_missense,trait){
  coding_sig$STAARB[coding_sig$STAARB == 0] <- 10^(-250)
  coding_minp <- min(coding_sig$STAARB)
  min_y <- ceiling(-log10(coding_minp)) + 1
  cex_point <- 1
  
  
  observed <- sort(coding_sig_plof$STAARB)
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  observed <- sort(coding_sig_plofds$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_missense$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### disruptive_missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_disruptive_missense$STAARB)
  
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(coding_sig_synonymous$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),ncol=1,bty="o",box.lwd=1,pch=0:4,cex=1.3,text.font=2)
  
}

print("WES Coding STAARB")

png(paste0("UKB_Imputed_Continuous_RV_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), 2, 6, byrow = TRUE))

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]
  
  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}

dev.off()

png(paste0("UKB_Imputed_Binary_RV_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2, 6, byrow = TRUE))

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/GeneCentricCoding/",trait,"_coding_sig.csv"))
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]
  
  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}
dev.off()


rm(list = ls())

coding_qqplot_STAARB <- function(coding_sig,coding_sig_plof,coding_sig_plofds,coding_sig_missense,coding_sig_synonymous,coding_sig_disruptive_missense,trait){
  coding_sig$STAARB[coding_sig$STAARB == 0] <- 10^(-250)
  coding_minp <- min(coding_sig$STAARB)
  min_y <- ceiling(-log10(coding_minp)) + 1
  cex_point <- 1
  
  
  observed <- sort(coding_sig_plof$STAARB)
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  observed <- sort(coding_sig_plofds$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_missense$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### disruptive_missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_disruptive_missense$STAARB)
  
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(coding_sig_synonymous$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),ncol=1,bty="o",box.lwd=1,pch=0:4,cex=1.3,text.font=2)
  
}

png(paste0("UKB_WGS_Binary_RV_Coding_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2, 6, byrow = TRUE))

print("WGS Coding STAARB Binary")

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/coding_sig.csv"))
  coding_sig <- coding_sig[coding_sig$Trait == trait,]
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]

  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}
dev.off()

print("WGS Coding STAARB Continuous")

png(paste0("UKB_WGS_Continuous_RV_Coding_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), 2, 6, byrow = TRUE))

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/coding_sig.csv"))
  coding_sig <- coding_sig[coding_sig$Trait == trait,]
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]
  
  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}
dev.off()


rm(list = ls())

noncoding_qqplot_STAARB <- function(noncoding_sig,noncoding_sig_upstream,noncoding_sig_downstream,noncoding_sig_ncRNA,noncoding_sig_UTR,noncoding_sig_promoter_CAGE,noncoding_sig_promoter_DHS,noncoding_sig_enhancer_CAGE,noncoding_sig_enhancer_DHS,trait){
  noncoding_minp <- min(noncoding_sig$STAARB)
  min_y <- ceiling(-log10(noncoding_minp)) + 1
  cex_point <- 1
  
  
  observed <- sort(noncoding_sig_upstream$STAARB)
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  
  observed <- sort(noncoding_sig_downstream$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### missense
  ## remove unconverged p-values
  observed <- sort(noncoding_sig_ncRNA$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### disruptive_missense
  ## remove unconverged p-values
  observed <- sort(noncoding_sig_UTR$STAARB)
  
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(noncoding_sig_promoter_CAGE$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(noncoding_sig_promoter_DHS$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=5, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(noncoding_sig_enhancer_CAGE$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=6, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(noncoding_sig_enhancer_DHS$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=7, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  legend("topleft",legend=c("upstream","downstream","ncRNA","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),ncol=1,bty="o",box.lwd=1,pch=0:7,cex=1,text.font=2)
  
}

print("WGS Noncoding STAARB Binary")

png(paste0("UKB_WGS_Binary_RV_Noncoding_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), 2, 6, byrow = TRUE))

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/noncoding_sig.csv"))
  noncoding_sig <- noncoding_sig[noncoding_sig$Trait == trait,]
  
  noncoding_sig_UTR <- noncoding_sig[noncoding_sig$Category == "UTR",]
  noncoding_sig_promoter_CAGE <- noncoding_sig[noncoding_sig$Category == "promoter_CAGE",]
  noncoding_sig_promoter_DHS <- noncoding_sig[noncoding_sig$Category == "promoter_DHS",]
  noncoding_sig_enhancer_CAGE <- noncoding_sig[noncoding_sig$Category == "enhancer_CAGE",]
  noncoding_sig_enhancer_DHS <- noncoding_sig[noncoding_sig$Category == "enhancer_DHS",]
  noncoding_sig_upstream <- noncoding_sig[noncoding_sig$Category == "upstream",]
  noncoding_sig_downstream <- noncoding_sig[noncoding_sig$Category == "downstream",]
  noncoding_sig_ncRNA <- noncoding_sig[noncoding_sig$Category == "ncRNA",]
  
  noncoding_qqplot_STAARB(noncoding_sig = noncoding_sig, noncoding_sig_ncRNA = noncoding_sig_ncRNA, noncoding_sig_downstream = noncoding_sig_downstream, noncoding_sig_upstream = noncoding_sig_upstream, noncoding_sig_UTR = noncoding_sig_UTR,noncoding_sig_promoter_CAGE = noncoding_sig_promoter_CAGE, noncoding_sig_promoter_DHS = noncoding_sig_promoter_DHS,noncoding_sig_enhancer_CAGE = noncoding_sig_enhancer_CAGE,noncoding_sig_enhancer_DHS = noncoding_sig_enhancer_DHS,trait = trait) 
}
dev.off()

print("WGS Noncoding STAARB Continuous")

png(paste0("UKB_WGS_Continuous_RV_Noncoding_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), 2, 6, byrow = TRUE))

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  noncoding_sig <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/noncoding_sig.csv"))
  noncoding_sig <- noncoding_sig[noncoding_sig$Trait == trait,]
  
  noncoding_sig_UTR <- noncoding_sig[noncoding_sig$Category == "UTR",]
  noncoding_sig_promoter_CAGE <- noncoding_sig[noncoding_sig$Category == "promoter_CAGE",]
  noncoding_sig_promoter_DHS <- noncoding_sig[noncoding_sig$Category == "promoter_DHS",]
  noncoding_sig_enhancer_CAGE <- noncoding_sig[noncoding_sig$Category == "enhancer_CAGE",]
  noncoding_sig_enhancer_DHS <- noncoding_sig[noncoding_sig$Category == "enhancer_DHS",]
  noncoding_sig_upstream <- noncoding_sig[noncoding_sig$Category == "upstream",]
  noncoding_sig_downstream <- noncoding_sig[noncoding_sig$Category == "downstream",]
  noncoding_sig_ncRNA <- noncoding_sig[noncoding_sig$Category == "ncRNA",]
  
  noncoding_qqplot_STAARB(noncoding_sig = noncoding_sig, noncoding_sig_ncRNA = noncoding_sig_ncRNA, noncoding_sig_downstream = noncoding_sig_downstream, noncoding_sig_upstream = noncoding_sig_upstream, noncoding_sig_UTR = noncoding_sig_UTR,noncoding_sig_promoter_CAGE = noncoding_sig_promoter_CAGE, noncoding_sig_promoter_DHS = noncoding_sig_promoter_DHS,noncoding_sig_enhancer_CAGE = noncoding_sig_enhancer_CAGE,noncoding_sig_enhancer_DHS = noncoding_sig_enhancer_DHS,trait = trait) 
}
dev.off()


rm(list = ls())

coding_qqplot_STAARB <- function(coding_sig,coding_sig_plof,coding_sig_plofds,coding_sig_missense,coding_sig_synonymous,coding_sig_disruptive_missense,trait){
  coding_sig$STAARB[coding_sig$STAARB == 0] <- 10^(-250)
  coding_minp <- min(coding_sig$STAARB)
  min_y <- ceiling(-log10(coding_minp)) + 1
  cex_point <- 1
  
  
  observed <- sort(coding_sig_plof$STAARB)
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  observed <- sort(coding_sig_plofds$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_missense$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### disruptive_missense
  ## remove unconverged p-values
  observed <- sort(coding_sig_disruptive_missense$STAARB)
  
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  ### synonymous
  ## remove unconverged p-values
  observed <- sort(coding_sig_synonymous$STAARB)
  
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  par(new=T)
  # par(mar=c(5,6,4,4))
  plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab=2,cex.lab=1,cex.axis=1,font.axis=2)
  title(main = trait,cex.main = 2)
  
  abline(0, 1, col="red",lwd=1)
  
  legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),ncol=1,bty="o",box.lwd=1,pch=0:4,cex=1.3,text.font=2)
  
}

print("AoU Coding STAARB Binary")
png(paste0("AoU_WES_Continuous_RV_QQplot.png"), width=15, height=9,units = "in",res = 300)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), 2, 6, byrow = TRUE))

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  coding_sig <- read.csv(paste0("/data/williamsjacr/AoU_Results/",trait,"_coding_sig.csv"))
  
  coding_sig_plof <- coding_sig[coding_sig$Category == "plof",]
  coding_sig_plofds <- coding_sig[coding_sig$Category == "plof_ds",]
  coding_sig_missense <- coding_sig[coding_sig$Category == "missense",]
  coding_sig_synonymous <- coding_sig[coding_sig$Category == "synonymous",]
  coding_sig_disruptive_missense <- coding_sig[coding_sig$Category == "disruptive_missense",]
  
  coding_qqplot_STAARB(coding_sig = coding_sig, coding_sig_plof = coding_sig_plof,coding_sig_plofds = coding_sig_plofds, coding_sig_missense = coding_sig_missense,coding_sig_synonymous = coding_sig_synonymous,coding_sig_disruptive_missense = coding_sig_disruptive_missense,trait = trait) 
}
dev.off()
