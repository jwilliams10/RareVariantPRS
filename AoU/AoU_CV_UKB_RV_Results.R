# =============================================================================
# [RICE-ANNOTATION] AoU/AoU_CV_UKB_RV_Results.R
# Purpose: Aggregates results across traits/ancestries and generates summary outputs.
#
# Paper linkage:
#   - Manuscript: Results -> All of Us Results (Fig. 7) and AoU-trained PRS evaluated on UKB (Fig. 8).
#   - Supplementary Figures: Supp. Fig. 17–21 (AoU association diagnostics, PRS performance, and portability).
#   - Supplementary Data: key AoU cohort summaries and diagnostics (e.g., sample sizes/variant counts/GC).
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
rm(list = ls())
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)

full_results <- NULL

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV (UKB)"
  full_results <- rbind(full_results,RICE_CV_Results)
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV (UKB)"
  full_results <- rbind(full_results,RICE_RV_Results)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CV_UKB_RV/",trait,"Best_Betas_RICECV.csv"))
  tmp$Method <- "RICE-CV (AoU)"
  full_results <- rbind(full_results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CV_UKB_RV/",trait,"Best_Betas_RICERV.csv"))
  tmp$Method <- "RICE-RV (AoU)"
  full_results <- rbind(full_results,tmp)
}

full_results <- full_results[full_results$ancestry %in% c("AFR","AMR","EUR","SAS"),]

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

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
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0



ylim <- max(c(full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("RICE-CV (UKB)","RICE-CV (AoU)","RICE-RV (UKB)","RICE-RV (AoU)"))

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) +
  ylab("Beta of PRS per SD") +
  ylim(0,ylim) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="Method + Dataset"))

ggsave(paste0("AoU_CrossPlatform_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)