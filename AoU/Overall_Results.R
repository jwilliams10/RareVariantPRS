rm(list = ls())

full_results <- read.csv("~/Desktop/RareVariantPRS_Results/AoU_Results.csv")

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
#full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CTSLEB","JointPRS","PROSPER","RICE-CV","RICE-RV"),]
full_results$Method[full_results$Method == "CTSLEB"] <- "CT-SLEB"


full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EAS","EUR","MID","SAS"))

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

pdf(paste0("Desktop/RareVariantPRS_Results/Figures/AoU_WES_Continuous_Raw_Beta.pdf"), width=10, height=6.18047)

ggplot(full_results) +
  geom_bar(aes(x=Method1, y=abs(beta_raw),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) +
  ggtitle("AoU WES Raw PRS Results") +
  ylab("Beta of PRS per SD") +
  ylim(0,0.7) +
  theme_Publication() +
  scale_fill_Publication()

dev.off()

# ggplot(full_results) +
#   geom_bar(aes(x=Method, y=abs(beta_raw),fill=Method), stat="identity", alpha=0.7) +
#   facet_grid(vars(trait), vars(ancestry)) + 
#   ggtitle("AoU WGS Raw PRS Results") + 
#   ylab("Beta of PRS per SD") + 
#   ylim(0,0.6) +
#   theme_Publication() + 
#   scale_fill_Publication()

pdf(paste0("Desktop/RareVariantPRS_Results/Figures/AoU_WES_Continuous_Adjusted_Beta.pdf"), width=10, height=6.18047)

ggplot(full_results) +
  geom_bar(aes(x=Method1, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +
  facet_grid(vars(trait), vars(ancestry)) +
  ggtitle("AoU WES Ancestry Adjusted PRS Results") +
  ylab("Beta of PRS per SD") +
  ylim(0,0.7) +
  theme_Publication() +
  scale_fill_Publication()

dev.off()

# ggplot(full_results) +
#   geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
#   # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
#   facet_grid(vars(trait), vars(ancestry)) + 
#   ggtitle("AoU WGS Ancestry Adjusted PRS Results") + 
#   ylab("Beta of PRS per SD") + 
#   ylim(0,0.6) +
#   theme_Publication() + 
#   scale_fill_Publication()