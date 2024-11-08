rm(list = ls())

# full_results <- NULL
# 
# for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
#   CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"Best_Betas.csv"))
#   CT_Results$Method <- "CT"
#   LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"Best_Betas.csv"))
#   LDPred2_Results$Method <- "LDPred"
#   LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"Best_Betas.csv"))
#   LASSOSUM2_Results$Method <- "LASSOSum"
#   CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"Best_Betas.csv"))
#   CV_Results$Method <- "CV_SL"
#   CV_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Common_plus_RareVariants/",trait,"Best_Betas.csv"))
#   full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,CV_Results,CV_RV_Results))
# }
# 
# full_results$beta_adjusted[full_results$beta_adjusted < 0 & full_results$Method %in% c("LDPred","LASSOSum")] <- -1*full_results$beta_adjusted[full_results$beta_adjusted < 0 & full_results$Method %in% c("LDPred","LASSOSum")]
# full_results$beta_raw[full_results$beta_raw < 0 & full_results$Method %in% c("LDPred","LASSOSum")] <- -1*full_results$beta_raw[full_results$beta_raw < 0 & full_results$Method %in% c("LDPred","LASSOSum")]
# 
# full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
# full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Continuous.csv")

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
full_results <- full_results[full_results$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results <- full_results[full_results$Method %in% c("CT","LASSOSum","LDPred","CV","RV"),]
full_results$Method[full_results$Method == "CV"] <- "RICE-CV" 
full_results$Method[full_results$Method == "RV"] <- "RICE-RV" 
full_results$Method[full_results$Method == "LDPred"] <- "LDpred2"
full_results$Method[full_results$Method == "LASSOSum"] <- "Lassosum2"


full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower = full_results$beta_raw - 1.96*full_results$se_raw, upper = full_results$beta_raw + 1.96*full_results$se_raw,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower = full_results$beta_adjusted - 1.96*full_results$se_adjusted, upper = full_results$beta_adjusted + 1.96*full_results$se_adjusted,method = full_results$Method,Standardization = "Using PCs 1-5"))

g1 <- ggplot(data=full_results_stacked[full_results_stacked$trait %in% c("BMI","Height"),], aes(x=method, y=beta, ymin=lower, ymax=upper,color = Standardization)) +
  geom_pointrange(position=position_dodge(width=.25),size = 0.2) + 
  facet_grid(vars(trait), vars(ancestry), scales="free") + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Beta of PRS per SD") +
  theme_Publication() + 
  scale_fill_Publication()
ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_Continuous_","A","_Raw_vs_AncestryAdjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)

g1 <- ggplot(data=full_results_stacked[full_results_stacked$trait %in% c("HDL","LDL"),], aes(x=method, y=beta, ymin=lower, ymax=upper,color = Standardization)) +
  geom_pointrange(position=position_dodge(width=.25),size = 0.2) + 
  facet_grid(vars(trait), vars(ancestry), scales="free") + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Beta of PRS per SD") +
  theme_Publication() + 
  scale_fill_Publication()
ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_Continuous_","B","_Raw_vs_AncestryAdjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)


g1 <- ggplot(data=full_results_stacked[full_results_stacked$trait %in% c("log(TG)","TC"),], aes(x=method, y=beta, ymin=lower, ymax=upper,color = Standardization)) +
  geom_pointrange(position=position_dodge(width=.25),size = 0.2) + 
  facet_grid(vars(trait), vars(ancestry), scales="free") + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Beta of PRS per SD") +
  theme_Publication() + 
  scale_fill_Publication()
ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_Continuous_","C","_Raw_vs_AncestryAdjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)


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

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method1, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("UKB WES PRS Results for Six Continuous Traits") + 
  ylab("Beta of PRS per SD") + 
  ylim(0,0.7) +
  theme_Publication() + 
  scale_fill_Publication()

ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_Continuous_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)