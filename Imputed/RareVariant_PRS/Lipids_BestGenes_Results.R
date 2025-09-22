rm(list = ls())

full_results <- NULL

for(trait in c("TC","HDL","LDL","logTG")){
  # tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_Lipids_Best_Best_Betas.csv"))
  # tmp$Method <- "High-Penetrance Genes (Functional Categories)"
  # tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  # full_results <- rbind(full_results,tmp) 
  
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_Lipids_Gene_Best_Betas.csv"))
  tmp$Method <- "High-Penetrance Genes"
  tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  full_results <- rbind(full_results,tmp) 
  
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"Best_Betas.csv"))
  tmp$Method <- "RICE-RV"
  tmp <- tmp[,c("trait","ancestry","Method","beta_adjusted")]
  full_results <- rbind(full_results,tmp) 
}

library(ggplot2)
library(boot)
full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"

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

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0

ylim <- max(full_results$beta_adjusted) + 0.05

full_results <- full_results[full_results$trait != "CAD",]

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("UKB Imputed + WES RICE-RV Results for Four Traits") + 
  ylab("Beta of PRS per SD") + 
  ylim(0,ylim) +
  theme_Publication() + 
  scale_fill_Publication()

ggsave(paste0("UKB_Imputed_BestGenes.png"),g2,width=10, height=6.18047,dpi = 300)
