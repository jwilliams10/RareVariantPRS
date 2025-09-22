rm(list = ls())

library(dplyr)
library(ggplot2)
library(cowplot)

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            strip.text.x = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.line = element_line(colour="black",size=2),
            axis.text.y = element_text(size = 12),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.text=element_text(size=12),
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

trait <- "Height"

for(trait in c("BMI","HDL","LDL","logTG","TC","Height")){
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
  pheno_validation$y_validation <- NA
  pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
  
  CV_RV_PRS_raw <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for(i in c("RV_PRS","CV_PRS")){
    tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      CV_RV_PRS_adjusted[,i] <- 0
    }else{
      CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
    }
  }
  
  
  CV_RV_PRS_adjusted$Common_Bin <- 9
  CV_RV_PRS_adjusted$Rare_Bin <- 9
  
  Common_quants <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(0,.1,.2,.3,.4,.6,.7,.8,.9))
  
  for(i in 1:8){
    CV_RV_PRS_adjusted$Common_Bin[CV_RV_PRS_adjusted$CV_PRS < unname(Common_quants[i + 1]) & CV_RV_PRS_adjusted$CV_PRS >= unname(Common_quants[i])] <- i
  }
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_EUR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_EUR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_EUR$y_validation <- scale(CV_RV_PRS_adjusted_EUR$y_validation)
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted_EUR[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin <- factor(CV_RV_PRS_adjusted_EUR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_AMR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AMR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_AMR$y_validation <- scale(CV_RV_PRS_adjusted_AMR$y_validation)
  
  CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted_AMR[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_AMR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AMR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_AFR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AFR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_AFR$y_validation <- scale(CV_RV_PRS_adjusted_AFR$y_validation)
  
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted_AFR[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_AFR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AFR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_SAS$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_SAS$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_SAS$y_validation <- scale(CV_RV_PRS_adjusted_SAS$y_validation)
  
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted_SAS[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_SAS$Rare_Bin <- factor(CV_RV_PRS_adjusted_SAS$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  CV_RV_PRS_adjusted_EUR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_EUR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_EUR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_EUR,mean)
  
  colnames(CV_RV_PRS_adjusted_EUR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_EUR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_EUR <- inner_join(CV_RV_PRS_adjusted_EUR,CV_RV_PRS_adjusted_EUR_se)
  
  CV_RV_PRS_adjusted_AMR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AMR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_AMR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AMR,mean)
  
  colnames(CV_RV_PRS_adjusted_AMR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_AMR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_AMR <- inner_join(CV_RV_PRS_adjusted_AMR,CV_RV_PRS_adjusted_AMR_se)
  
  CV_RV_PRS_adjusted_AFR_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AFR,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_AFR <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_AFR,mean)
  
  colnames(CV_RV_PRS_adjusted_AFR) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_AFR_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_AFR <- inner_join(CV_RV_PRS_adjusted_AFR,CV_RV_PRS_adjusted_AFR_se)
  
  CV_RV_PRS_adjusted_SAS_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_SAS,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted_SAS <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted_SAS,mean)
  
  colnames(CV_RV_PRS_adjusted_SAS) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_SAS_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted_SAS <- inner_join(CV_RV_PRS_adjusted_SAS,CV_RV_PRS_adjusted_SAS_se)
  
  
  colnames(CV_RV_PRS_adjusted_EUR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_AMR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_AFR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  colnames(CV_RV_PRS_adjusted_SAS) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  
  ymin_EUR <- round(min(c(CV_RV_PRS_adjusted_EUR$Mean - CV_RV_PRS_adjusted_EUR$SE)) - 0.05,2)
  ymax_EUR <- round(max(c(CV_RV_PRS_adjusted_EUR$Mean + CV_RV_PRS_adjusted_EUR$SE)) + 0.05,2)
  
  ymin_AMR <- round(min(c(CV_RV_PRS_adjusted_AMR$Mean - CV_RV_PRS_adjusted_AMR$SE)) - 0.05,2)
  ymax_AMR <- round(max(c(CV_RV_PRS_adjusted_AMR$Mean + CV_RV_PRS_adjusted_AMR$SE)) + 0.05,2)
  
  ymin_AFR <- round(min(c(CV_RV_PRS_adjusted_AFR$Mean - CV_RV_PRS_adjusted_AFR$SE)) - 0.05,2)
  ymax_AFR <- round(max(c(CV_RV_PRS_adjusted_AFR$Mean + CV_RV_PRS_adjusted_AFR$SE)) + 0.05,2)
  
  ymin_SAS <- round(min(c(CV_RV_PRS_adjusted_SAS$Mean - CV_RV_PRS_adjusted_SAS$SE)) - 0.05,2)
  ymax_SAS <- round(max(c(CV_RV_PRS_adjusted_SAS$Mean + CV_RV_PRS_adjusted_SAS$SE)) + 0.05,2)
  
  plot1 <- ggplot(data=CV_RV_PRS_adjusted_EUR, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0("Standardized ",ifelse(trait == "logTG","log(TG)",trait))) + ylim(c(ymin_EUR,ymax_EUR)) +
    scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%")) + labs(x = "RICE-CV Quantiles (Common Variants)")
  plot2 <- ggplot(data=CV_RV_PRS_adjusted_AMR, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0("Standardized ",ifelse(trait == "logTG","log(TG)",trait))) + ylim(c(ymin_AMR,ymax_AMR)) +
    scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%")) + labs(x = "RICE-CV Quantiles (Common Variants)")
  plot3 <- ggplot(data=CV_RV_PRS_adjusted_AFR, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0("Standardized ",ifelse(trait == "logTG","log(TG)",trait))) + ylim(c(ymin_AFR,ymax_AFR)) +
    scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%")) + labs(x = "RICE-CV Quantiles (Common Variants)")
  plot4 <- ggplot(data=CV_RV_PRS_adjusted_SAS, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0("Standardized ",ifelse(trait == "logTG","log(TG)",trait))) + ylim(c(ymin_SAS,ymax_SAS)) +
    scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%")) + labs(x = "RICE-CV Quantiles (Common Variants)")
  
  
  prow <- plot_grid(
    plot1 + theme(legend.position="none"),
    plot2 + theme(legend.position="none"),
    plot3 + theme(legend.position="none"),
    plot4 + theme(legend.position="none"),
    align = 'vh',
    labels = c("EUR","AMR","AFR","SAS"),
    hjust = -1,
    ncol = 2,label_size = 16
  )
  
  legend_b <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]
  
  pdf(paste0(trait,"_Imputed_Fig5.pdf"), width=12, height=7.416564)
  
  print(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
  
  dev.off()
}


















rm(list = ls())

library(dplyr)
library(ggplot2)
library(cowplot)

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            strip.text.x = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.line = element_line(colour="black",size=2),
            axis.text.y = element_text(size = 12),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.text=element_text(size=12),
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

trait <- "T2D"

for(trait in c("Asthma","T2D","CAD","Breast","Prostate")){
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  CV_RV_PRS_raw <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for(i in c("RV_PRS","CV_PRS")){
    tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    R <- mod$residuals
    tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
    mod <- lm(y~.,data = tmp)
    y_hat <- predict(mod,tmp)
    if(sum(y_hat < 0) > 0){
      mod <- lm(y~1,data = tmp)
      y_hat <- predict(mod,tmp)
    }
    if(sum(sqrt(y_hat)) == 0){
      CV_RV_PRS_adjusted[,i] <- 0
    }else{
      CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
    }
  }
  
  
  CV_RV_PRS_adjusted$Common_Bin <- 5
  CV_RV_PRS_adjusted$Rare_Bin <- 9
  
  Common_quants <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(0,.2,.4,.6,.8))
  
  for(i in 1:4){
    CV_RV_PRS_adjusted$Common_Bin[CV_RV_PRS_adjusted$CV_PRS < unname(Common_quants[i + 1]) & CV_RV_PRS_adjusted$CV_PRS >= unname(Common_quants[i])] <- i
  }
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  # CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
  # CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  # CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  # 
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted_EUR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_EUR$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted_EUR$Rare_Bin[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted_EUR[CV_RV_PRS_adjusted_EUR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted_EUR$Rare_Bin <- factor(CV_RV_PRS_adjusted_EUR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  # Rare_quants <- quantile(CV_RV_PRS_adjusted_AMR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  # for(i in 1:8){
  #   CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AMR$RV_PRS >= unname(Rare_quants[i])] <- i
  # }
  # 
  # CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin == 1] <- "Below 5%"
  # CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  # CV_RV_PRS_adjusted_AMR$Rare_Bin[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("9")] <- "Above 95%"
  # 
  # CV_RV_PRS_adjusted_AMR <- CV_RV_PRS_adjusted_AMR[CV_RV_PRS_adjusted_AMR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  # 
  # CV_RV_PRS_adjusted_AMR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AMR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  # 
  # Rare_quants <- quantile(CV_RV_PRS_adjusted_AFR$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  # for(i in 1:8){
  #   CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_AFR$RV_PRS >= unname(Rare_quants[i])] <- i
  # }
  # 
  # CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin == 1] <- "Below 5%"
  # CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  # CV_RV_PRS_adjusted_AFR$Rare_Bin[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("9")] <- "Above 95%"
  # 
  # CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted_AFR[CV_RV_PRS_adjusted_AFR$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  # 
  # CV_RV_PRS_adjusted_AFR$Rare_Bin <- factor(CV_RV_PRS_adjusted_AFR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  # 
  # Rare_quants <- quantile(CV_RV_PRS_adjusted_SAS$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  # for(i in 1:8){
  #   CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted_SAS$RV_PRS >= unname(Rare_quants[i])] <- i
  # }
  # 
  # CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin == 1] <- "Below 5%"
  # CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  # CV_RV_PRS_adjusted_SAS$Rare_Bin[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("9")] <- "Above 95%"
  # 
  # CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted_SAS[CV_RV_PRS_adjusted_SAS$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  # 
  # CV_RV_PRS_adjusted_SAS$Rare_Bin <- factor(CV_RV_PRS_adjusted_SAS$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  log_odds_ratio_EUR <- NULL
  for(rv_bin in c("Below 5%","30% - 70%","Above 95%")){
    for(cv_bin in 1:5){
      
      if(rv_bin == "30% - 70%" & cv_bin == 3){
        log_odds_ratio <- log(1)
        se_odds <- 0
      }
      
      tmp <- data.frame(Group = c(rep("Control",sum(CV_RV_PRS_adjusted_EUR$Common_Bin == 3 & CV_RV_PRS_adjusted_EUR$Rare_Bin == "30% - 70%")),rep("Comparison",sum(CV_RV_PRS_adjusted_EUR$Common_Bin == cv_bin & CV_RV_PRS_adjusted_EUR$Rare_Bin == rv_bin))),Case_Control = c(CV_RV_PRS_adjusted_EUR[CV_RV_PRS_adjusted_EUR$Common_Bin == 3 & CV_RV_PRS_adjusted_EUR$Rare_Bin == "30% - 70%",trait],CV_RV_PRS_adjusted_EUR[CV_RV_PRS_adjusted_EUR$Common_Bin == cv_bin & CV_RV_PRS_adjusted_EUR$Rare_Bin == rv_bin,trait]))
      
      log_odds_ratio <- log((table(tmp)[1,2]/table(tmp)[1,1])/(table(tmp)[2,2]/table(tmp)[2,1]))
      se_log_odds_ratio <- sqrt(1/table(tmp)[1,2] + 1/table(tmp)[1,1] + 1/table(tmp)[2,2] + 1/table(tmp)[2,1])
      
      log_odds_ratio_EUR <- rbind(log_odds_ratio_EUR,data.frame(Common_Bin = cv_bin,Rare_Bin = rv_bin,log_odds_ratio = log_odds_ratio,se_log_odds_ratio = se_log_odds_ratio))
    }
  }
  
  # log_odds_ratio_AFR <- NULL
  # for(rv_bin in c("Below 5%","30% - 70%","Above 95%")){
  #   for(cv_bin in 1:9){
  #     
  #     if(rv_bin == "30% - 70%" & cv_bin == 5){
  #       log_odds_ratio <- log(1)
  #       se_odds <- 0
  #     }
  #     
  #     tmp <- data.frame(Group = c(rep("Control",sum(CV_RV_PRS_adjusted_AFR$Common_Bin == 5 & CV_RV_PRS_adjusted_AFR$Rare_Bin == "30% - 70%")),rep("Comparison",sum(CV_RV_PRS_adjusted_AFR$Common_Bin == cv_bin & CV_RV_PRS_adjusted_AFR$Rare_Bin == rv_bin))),Case_Control = c(CV_RV_PRS_adjusted_AFR$y_validation[CV_RV_PRS_adjusted_AFR$Common_Bin == 5 & CV_RV_PRS_adjusted_AFR$Rare_Bin == "30% - 70%"],CV_RV_PRS_adjusted_AFR$y_validation[CV_RV_PRS_adjusted_AFR$Common_Bin == cv_bin & CV_RV_PRS_adjusted_AFR$Rare_Bin == rv_bin]))
  #     
  #     log_odds_ratio <- log((table(tmp)[1,2]/table(tmp)[1,1])/(table(tmp)[2,2]/table(tmp)[2,1]))
  #     se_log_odds_ratio <- sqrt(1/table(tmp)[1,2] + 1/table(tmp)[1,1] + 1/table(tmp)[2,2] + 1/table(tmp)[2,1])
  #     
  #     log_odds_ratio_AFR <- rbind(log_odds_ratio_AFR,data.frame(Common_Bin = cv_bin,Rare_Bin = rv_bin,log_odds_ratio = log_odds_ratio,se_log_odds_ratio = se_log_odds_ratio))
  #   }
  # }
  # 
  # log_odds_ratio_SAS <- NULL
  # for(rv_bin in c("Below 5%","30% - 70%","Above 95%")){
  #   for(cv_bin in 1:9){
  #     
  #     if(rv_bin == "30% - 70%" & cv_bin == 5){
  #       log_odds_ratio <- log(1)
  #       se_odds <- 0
  #     }
  #     
  #     tmp <- data.frame(Group = c(rep("Control",sum(CV_RV_PRS_adjusted_SAS$Common_Bin == 5 & CV_RV_PRS_adjusted_SAS$Rare_Bin == "30% - 70%")),rep("Comparison",sum(CV_RV_PRS_adjusted_SAS$Common_Bin == cv_bin & CV_RV_PRS_adjusted_SAS$Rare_Bin == rv_bin))),Case_Control = c(CV_RV_PRS_adjusted_SAS$y_validation[CV_RV_PRS_adjusted_SAS$Common_Bin == 5 & CV_RV_PRS_adjusted_SAS$Rare_Bin == "30% - 70%"],CV_RV_PRS_adjusted_SAS$y_validation[CV_RV_PRS_adjusted_SAS$Common_Bin == cv_bin & CV_RV_PRS_adjusted_SAS$Rare_Bin == rv_bin]))
  #     
  #     log_odds_ratio <- log((table(tmp)[1,2]/table(tmp)[1,1])/(table(tmp)[2,2]/table(tmp)[2,1]))
  #     se_log_odds_ratio <- sqrt(1/table(tmp)[1,2] + 1/table(tmp)[1,1] + 1/table(tmp)[2,2] + 1/table(tmp)[2,1])
  #     
  #     log_odds_ratio_SAS <- rbind(log_odds_ratio_SAS,data.frame(Common_Bin = cv_bin,Rare_Bin = rv_bin,log_odds_ratio = log_odds_ratio,se_log_odds_ratio = se_log_odds_ratio))
  #   }
  # }
  # 
  # log_odds_ratio_AMR <- NULL
  # for(rv_bin in c("Below 5%","30% - 70%","Above 95%")){
  #   for(cv_bin in 1:9){
  #     
  #     if(rv_bin == "30% - 70%" & cv_bin == 5){
  #       log_odds_ratio <- log(1)
  #       se_odds <- 0
  #     }
  #     
  #     tmp <- data.frame(Group = c(rep("Control",sum(CV_RV_PRS_adjusted_AMR$Common_Bin == 5 & CV_RV_PRS_adjusted_AMR$Rare_Bin == "30% - 70%")),rep("Comparison",sum(CV_RV_PRS_adjusted_AMR$Common_Bin == cv_bin & CV_RV_PRS_adjusted_AMR$Rare_Bin == rv_bin))),Case_Control = c(CV_RV_PRS_adjusted_AMR$y_validation[CV_RV_PRS_adjusted_AMR$Common_Bin == 5 & CV_RV_PRS_adjusted_AMR$Rare_Bin == "30% - 70%"],CV_RV_PRS_adjusted_AMR$y_validation[CV_RV_PRS_adjusted_AMR$Common_Bin == cv_bin & CV_RV_PRS_adjusted_AMR$Rare_Bin == rv_bin]))
  #     
  #     print(table(tmp))
  #     
  #     log_odds_ratio <- log((table(tmp)[1,2]/table(tmp)[1,1])/(table(tmp)[2,2]/table(tmp)[2,1]))
  #     se_log_odds_ratio <- sqrt(1/table(tmp)[1,2] + 1/table(tmp)[1,1] + 1/table(tmp)[2,2] + 1/table(tmp)[2,1])
  #     
  #     log_odds_ratio_AMR <- rbind(log_odds_ratio_AMR,data.frame(Common_Bin = cv_bin,Rare_Bin = rv_bin,log_odds_ratio = log_odds_ratio,se_log_odds_ratio = se_log_odds_ratio))
  #   }
  # }
  # 
  # 
  log_odds_ratio_EUR$Rare_Bin <- factor(log_odds_ratio_EUR$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  colnames(log_odds_ratio_EUR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  # colnames(log_odds_ratio_AMR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  # colnames(log_odds_ratio_SAS) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  # colnames(log_odds_ratio_AFR) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  
  ymin <- round(min(exp(log_odds_ratio_EUR$Mean - log_odds_ratio_EUR$SE)) - 0.05,2)
  ymax <- round(max(exp(log_odds_ratio_EUR$Mean + log_odds_ratio_EUR$SE)) + 0.05,2)
  
  plot1 <- ggplot(data=log_odds_ratio_EUR, aes(x=Common_Bin, y=exp(Mean), color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=exp(Mean-SE), ymax=exp(Mean+SE))) + theme_Publication() + ylab(paste0(trait," Odds Ratio")) + ylim(c(ymin,ymax)) + 
    scale_x_continuous(breaks = c(1:5),labels = c("0-20%","20-40%","40-60%","60-80%","80-100%")) + labs(x = "RICE-CV Quantiles (Common Variants)")
  # plot2 <- ggplot(data=log_odds_ratio_AMR, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0(trait," Standardized")) + ylim(c(ymin,ymax)) + 
  #   scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%"))
  # plot3 <- ggplot(data=log_odds_ratio_AFR, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0(trait," Standardized")) + ylim(c(ymin,ymax)) + 
  #   scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%"))
  # plot4 <- ggplot(data=log_odds_ratio_SAS, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0(trait," Standardized")) + ylim(c(ymin,ymax)) + 
  #   scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%"))
  # 
  
  prow <- plot_grid(
    plot1 + theme(legend.position="none"),
    # plot2 + theme(legend.position="none"),
    # plot3 + theme(legend.position="none"),
    # plot4 + theme(legend.position="none"),
    align = 'vh',
    labels = c("EUR"),
    hjust = -1,
    ncol = 1,label_size = 16
  )
  
  legend_b <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]
  
  pdf(paste0(trait,"_Imputed_Fig5.pdf"), width=12, height=7.416564)
  
  print(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))
  
  dev.off()
}
