rm(list = ls())

full_results <- NULL

options(scipen = 0)

for(thresholds in c(1e-5,1e-4,1e-3,1e-2)){
  for(trait in c("BMI","LDL","HDL","logTG","TC","Height","Breast","Prostate","CAD","T2D","Asthma")){
    tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Threshold_Sensitivity_RICE_RV/",trait,"_",thresholds,"Best_Betas.csv"))
    tmp$Threshold <- thresholds
    tmp <- tmp[,c("Threshold","trait","ancestry","beta_adjusted")]
    full_results <- rbind(full_results,tmp) 
  }
}

options("scipen"=100, "digits"=4)
full_results$Threshold <- as.character(full_results$Threshold)
options(scipen = 0)

full_results$Threshold <- factor(full_results$Threshold,levels = c("0.00001","0.0001","0.001","0.01"))

library(ggplot2)
library(boot)
full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"

full_results_continuous <- full_results[full_results$trait %in% c("BMI","Height","HDL","LDL","log(TG)","TC"),]
full_results_binary <- full_results[full_results$trait %in% c("Asthma","CAD","T2D","Breast","Prostate"),]

full_results_continuous$trait <- factor(full_results_continuous$trait,levels = c("BMI","HDL","Height","LDL","log(TG)","TC"))
full_results_binary$trait <- factor(full_results_binary$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))

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

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#E23838","#973999","#009cdf")), ...)
  
}

full_results_continuous$beta_adjusted[full_results_continuous$beta_adjusted < 0] <- 0
full_results_binary$beta_adjusted[full_results_binary$beta_adjusted < 0] <- 0

ylim <- max(full_results_continuous$beta_adjusted) + 0.05

g2 <- ggplot(full_results_continuous) +
  geom_bar(aes(x=Threshold, y=abs(beta_adjusted),fill=Threshold), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("UKB Imputed + WES RICE-RV Results for Six Continuous Traits") + 
  ylab("Beta of PRS per SD") + 
  ylim(0,ylim) +
  theme_Publication() + 
  scale_fill_Publication()

ggsave(paste0("UKB_Imputed_Continuous_Thresholds.png"),g2,width=10, height=6.18047,dpi = 300)

ylim <- max(full_results_binary$beta_adjusted) + 0.05

g3 <- ggplot(full_results_binary) +
  geom_bar(aes(x=Threshold, y=abs(beta_adjusted),fill=Threshold), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("UKB Imputed + WES RICE-RV Results for Five Binary Traits") + 
  ylab("Beta of PRS per SD") + 
  ylim(0,ylim) +
  theme_Publication() + 
  scale_fill_Publication()

ggsave(paste0("UKB_Imputed_Binary_Thresholds.png"),g3,width=10, height=6.18047,dpi = 300)