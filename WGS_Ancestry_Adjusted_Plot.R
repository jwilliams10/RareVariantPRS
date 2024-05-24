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
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_blank(),
            # axis.text.x = element_blank(), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
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

plot_data <- NULL


for(trait in c("BMI","HDL","LDL","logTG","TC","Height")){
  RV_PRS <- read.csv(paste0("Desktop/RareVariantPRS_Results/Continuous_PRS_Validation/",trait,"_BestPRS.csv"))
  CV_PRS <- read.delim(paste0("Desktop/RareVariantPRS_Results/Continuous_PRS_Validation/",trait,"_Best_Validation_All.txt"))
  
  CV_RV_PRS <- inner_join(RV_PRS,CV_PRS)
  CV_RV_PRS_raw <- CV_RV_PRS
  CV_RV_PRS_adjusted <- CV_RV_PRS
  
  for(i in c("RV_PRS","prs")){
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
  
  load("Desktop/RareVariantPRS_Results/all_phenotypes.RData")
  
  CV_RV_PRS_raw_EUR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  CV_RV_PRS_raw_SAS <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  CV_RV_PRS_raw_MIX <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  CV_RV_PRS_raw_AFR <- CV_RV_PRS_raw[CV_RV_PRS_raw$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  
  CV_RV_PRS_adjusted_EUR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  CV_RV_PRS_adjusted_MIX <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  CV_RV_PRS_adjusted_AFR <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  CV_RV_PRS_adjusted_SAS <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  
  CV_RV_PRS_raw_EUR$Y <- scale(CV_RV_PRS_raw_EUR$Y)
  CV_RV_PRS_raw_MIX$Y <- scale(CV_RV_PRS_raw_MIX$Y)
  CV_RV_PRS_raw_AFR$Y <- scale(CV_RV_PRS_raw_AFR$Y)
  CV_RV_PRS_raw_SAS$Y <- scale(CV_RV_PRS_raw_SAS$Y)
  
  CV_RV_PRS_raw_EUR$Y <- scale(CV_RV_PRS_raw_EUR$Y)
  CV_RV_PRS_raw_MIX$Y <- scale(CV_RV_PRS_raw_MIX$Y)
  CV_RV_PRS_raw_AFR$Y <- scale(CV_RV_PRS_raw_AFR$Y)
  CV_RV_PRS_raw_SAS$Y <- scale(CV_RV_PRS_raw_SAS$Y)
  
  CV_RV_PRS_adjusted_EUR$Y <- scale(CV_RV_PRS_adjusted_EUR$Y)
  CV_RV_PRS_adjusted_MIX$Y <- scale(CV_RV_PRS_adjusted_MIX$Y)
  CV_RV_PRS_adjusted_AFR$Y <- scale(CV_RV_PRS_adjusted_AFR$Y)
  CV_RV_PRS_adjusted_SAS$Y <- scale(CV_RV_PRS_adjusted_SAS$Y)
  
  Raw <- data.frame(Y = c(CV_RV_PRS_raw_EUR$Y,CV_RV_PRS_raw_MIX$Y,CV_RV_PRS_raw_AFR$Y,CV_RV_PRS_raw_SAS$Y),
                    Ancestry = c(rep("EUR",length(CV_RV_PRS_raw_EUR$Y)),rep("MIX",length(CV_RV_PRS_raw_MIX$Y)),rep("AFR",length(CV_RV_PRS_raw_AFR$Y)),rep("SAS",length(CV_RV_PRS_raw_SAS$Y))),
                    Common_PRS = c(CV_RV_PRS_raw_EUR$prs,CV_RV_PRS_raw_MIX$prs,CV_RV_PRS_raw_AFR$prs,CV_RV_PRS_raw_SAS$prs),
                    Rare_PRS = c(CV_RV_PRS_raw_EUR$RV_PRS,CV_RV_PRS_raw_MIX$RV_PRS,CV_RV_PRS_raw_AFR$RV_PRS,CV_RV_PRS_raw_SAS$RV_PRS),
                    Raw = "Raw",Trait = trait)
  Adjusted <- data.frame(Y = c(CV_RV_PRS_adjusted_EUR$Y,CV_RV_PRS_adjusted_MIX$Y,CV_RV_PRS_adjusted_AFR$Y,CV_RV_PRS_adjusted_SAS$Y),
                         Ancestry = c(rep("EUR",length(CV_RV_PRS_adjusted_EUR$Y)),rep("MIX",length(CV_RV_PRS_adjusted_MIX$Y)),rep("AFR",length(CV_RV_PRS_adjusted_AFR$Y)),rep("SAS",length(CV_RV_PRS_adjusted_SAS$Y))),
                         Common_PRS = c(CV_RV_PRS_adjusted_EUR$prs,CV_RV_PRS_adjusted_MIX$prs,CV_RV_PRS_adjusted_AFR$prs,CV_RV_PRS_adjusted_SAS$prs),
                         Rare_PRS = c(CV_RV_PRS_adjusted_EUR$RV_PRS,CV_RV_PRS_adjusted_MIX$RV_PRS,CV_RV_PRS_adjusted_AFR$RV_PRS,CV_RV_PRS_adjusted_SAS$RV_PRS),
                         Raw = "Adjusted",Trait = trait)
  
  plot_data <- rbind(plot_data,rbind(Raw,Adjusted))
  
}

ggplot(plot_data,aes(x = Common_PRS,color = Ancestry)) +
  geom_density() +
  facet_grid(rows = vars(Trait),cols = vars(Raw)) + 
  ggtitle("WGS Ancestry Adjustment; RICE-CV") + 
  ylab("Density") + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(plot_data,aes(x = Rare_PRS,color = Ancestry)) +
  geom_density() +
  facet_grid(rows = vars(Trait),cols = vars(Raw), scales="free") + 
  ggtitle("WGS Ancestry Adjustment; RICE-RV") + 
  ylab("Density") + 
  theme_Publication() + 
  scale_fill_Publication()

plot_data <- plot_data[plot_data$Rare_PRS < quantile(plot_data$Rare_PRS,.99) & plot_data$Rare_PRS > quantile(plot_data$Rare_PRS,0.01),]
ggplot(plot_data,aes(x = Rare_PRS,color = Ancestry)) +
  geom_density() +
  facet_grid(rows = vars(Trait),cols = vars(Raw), scales="free") + 
  ggtitle("WGS Ancestry Adjustment; RICE-RV") + 
  ylab("Density") + 
  theme_Publication() + 
  scale_fill_Publication()
