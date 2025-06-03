rm(list = ls())

library(ggplot2)
library(ggpubr)
library(dplyr)

full_results_binary <- NULL

for(trait in c("Asthma","Breast","CAD","Prostate","T2D")){
  Coding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/",trait,"_Coding_Best_Betas.csv"))
  Coding_Results$Method <- "Coding"
  Noncoding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestRareVariantPRS/",trait,"_Noncoding_Best_Betas.csv"))
  Noncoding_Results$Method <- "Noncoding"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  full_results_binary <- rbind(full_results_binary,rbind(Coding_Results,Noncoding_Results,RICE_CV_Results,RICE_RV_Results))
}

# full_results_binary <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Binary.csv")

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
            axis.text.x = element_blank(), 
            axis.text.y = element_text(size = 10), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
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
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#E23838","#973999","#009cdf")), ...)
  
}


library(ggplot2)
library(boot)
full_results_binary <- full_results_binary[full_results_binary$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_binary <- full_results_binary[full_results_binary$Method %in% c("Coding","Noncoding","RICE-RV","RICE-CV"),]

full_results_binary$Method <- factor(full_results_binary$Method,levels = c("Coding","Noncoding","RICE-RV","RICE-CV"))
full_results_binary$trait <- factor(full_results_binary$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))
full_results_binary$ancestry <- factor(full_results_binary$ancestry,levels = c("AFR","AMR","EUR","SAS"))


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
            axis.text.x = element_blank(), 
            axis.text.y = element_text(size = 10), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
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

full_results_binary$beta_adjusted[full_results_binary$beta_adjusted < 0] <- 0
full_results_binary$beta_raw[full_results_binary$beta_raw < 0] <- 0

ylim <- max(full_results_binary$beta_adjusted) + 0.05

g2 <- ggplot(full_results_binary) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=AUC_low, ymax=AUC_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("UKB WGS PRS Results for Five Binary Traits") + 
  ylab("log(Odds Ratio) of PRS per SD") + 
  ylim(0,ylim) +
  theme_Publication() + 
  scale_fill_Publication()

# ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WGS_Binary_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)
ggsave(paste0("UKB_WGS_Binary_Coding_Noncoding_Beta.png"),g2,width=10, height=6.18047,dpi = 300)


rm(list = ls())

library(ggplot2)
library(ggpubr)
library(dplyr)

full_results_continuous <- NULL

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  Coding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"_Coding_Best_Betas.csv"))
  Coding_Results$Method <- "Coding"
  Noncoding_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"_Noncoding_Best_Betas.csv"))
  Noncoding_Results$Method <- "Noncoding"
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  full_results_continuous <- rbind(full_results_continuous,rbind(Coding_Results,Noncoding_Results,RICE_CV_Results,RICE_RV_Results))
}

# full_results_continuous <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Continuous.csv")

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
            axis.text.x = element_blank(), 
            axis.text.y = element_text(size = 10), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
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
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#E23838","#973999","#009cdf")), ...)
  
}


library(ggplot2)
library(boot)
full_results_continuous <- full_results_continuous[full_results_continuous$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_continuous <- full_results_continuous[full_results_continuous$Method %in% c("Coding","Noncoding","RICE-RV","RICE-CV"),]

full_results_continuous$Method <- factor(full_results_continuous$Method,levels = c("Coding","Noncoding","RICE-RV","RICE-CV"))
full_results_continuous$trait[full_results_continuous$trait == "logTG"] <- "log(TG)"
full_results_continuous$trait <- factor(full_results_continuous$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_continuous$ancestry <- factor(full_results_continuous$ancestry,levels = c("AFR","AMR","EUR","SAS"))


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
            axis.text.x = element_blank(), 
            axis.text.y = element_text(size = 10), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
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

full_results_continuous$beta_adjusted[full_results_continuous$beta_adjusted < 0] <- 0
full_results_continuous$beta_raw[full_results_continuous$beta_raw < 0] <- 0

ylim <- max(full_results_continuous$beta_adjusted) + 0.05

g2 <- ggplot(full_results_continuous) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=AUC_low, ymax=AUC_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("UKB WGS PRS Results for Six Continuous Traits") + 
  ylab("Beta of PRS per SD") + 
  ylim(0,ylim) +
  theme_Publication() + 
  scale_fill_Publication()

# ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WGS_Continuous_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)
ggsave(paste0("UKB_WGS_Continuous_Coding_Noncoding_Beta.png"),g2,width=10, height=6.18047,dpi = 300)