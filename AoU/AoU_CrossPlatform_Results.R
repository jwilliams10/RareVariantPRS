rm(list = ls())
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)

RICE_CV_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Results.csv")
RICE_CV_Results$Method <- "RICE-CV; Train (AoU) -> Validate (AoU)"
RICE_CV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Boot.csv")
RICE_CV_Boot_Results$Method <- "RICE-CV; Train (AoU) -> Validate (AoU)"
RICE_CV_Boot_Results <- RICE_CV_Boot_Results[,c("trait","beta_CV_adjusted_EUR_boot","beta_CV_adjusted_SAS_boot","beta_CV_adjusted_AMR_boot","beta_CV_adjusted_AFR_boot","Method")]
colnames(RICE_CV_Boot_Results) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")
RICE_RV_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Results.csv")
RICE_RV_Results$Method <- "RICE-RV; Train (AoU) -> Validate (AoU)"
RICE_RV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Boot.csv")
RICE_RV_Boot_Results$Method <- "RICE-RV; Train (AoU) -> Validate (AoU)"
RICE_RV_Boot_Results <- RICE_RV_Boot_Results[,c("trait","beta_RV_adjusted_EUR_boot","beta_RV_adjusted_SAS_boot","beta_RV_adjusted_AMR_boot","beta_RV_adjusted_AFR_boot","Method")]
colnames(RICE_RV_Boot_Results) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")

for(trait in c("BMI","LDL","HDL","logTG","TC","Height")){
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"Best_Betas_RICECV.csv"))
  tmp$Method <- "RICE-CV; Train (AoU) -> Validate (UKB)"
  RICE_CV_Results <- rbind(RICE_CV_Results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Bootstraps_RICECV.csv"))
  tmp$Method <- "RICE-CV; Train (AoU) -> Validate (UKB)"
  tmp <- tmp[,c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")]
  colnames(tmp) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")
  RICE_CV_Boot_Results <- rbind(RICE_CV_Boot_Results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"Best_Betas_RICERV.csv"))
  tmp$Method <- "RICE-RV; Train (AoU) -> Validate (UKB)"
  RICE_RV_Results <- rbind(RICE_RV_Results,tmp)
  tmp <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/AoU_CrossPlatform/",trait,"_Bootstraps_RICERV.csv"))
  tmp$Method <- "RICE-RV; Train (AoU) -> Validate (UKB)"
  tmp <- tmp[,c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")]
  colnames(tmp) <- c("trait","beta_adjusted_EUR_boot","beta_adjusted_SAS_boot","beta_adjusted_AMR_boot","beta_adjusted_AFR_boot","Method")
  RICE_RV_Boot_Results <- rbind(RICE_RV_Boot_Results,tmp)
}

full_results <- rbind(RICE_CV_Results,RICE_RV_Results)
full_results_Boot <- rbind(RICE_CV_Boot_Results,RICE_RV_Boot_Results)

full_results <- full_results[full_results$ancestry %in% c("AFR","AMR","EUR","SAS"),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("RICE-RV; Train (AoU) -> Validate (UKB)","RICE-CV; Train (AoU) -> Validate (UKB)","RICE-RV; Train (AoU) -> Validate (AoU)","RICE-CV; Train (AoU) -> Validate (AoU)"))
full_results$Method1[full_results$Method1 == "RICE-RV; Train (AoU) -> Validate (UKB)"] <- "RICE-CV; Train (AoU) -> Validate (UKB)"
full_results$Method1[full_results$Method1 == "RICE-RV; Train (AoU) -> Validate (AoU)"] <- "RICE-CV; Train (AoU) -> Validate (AoU)"
full_results$Method1 <- factor(full_results$Method1,levels = c("RICE-CV; Train (AoU) -> Validate (UKB)","RICE-CV; Train (AoU) -> Validate (AoU)"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EUR","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("AFR","AMR","EUR","SAS"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower,CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper,CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("AFR","AMR","EUR","SAS"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower,CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper,CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

full_results$Method_Dodge <- ifelse(str_detect(full_results$Method,"RICE-CV"),"RICE-CV","RICE-RV")

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

full_results$group1 <- "RICE-RV; Train (AoU) -> Validate (UKB)"
full_results$group2 <- "RICE-RV; Train (AoU) -> Validate (UKB)"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta1[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] > 0,"**",""))
full_results$position1 <- NA
full_results$position1[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] <- full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"] + 0.03

full_results$group3 <- "RICE-RV; Train (AoU) -> Validate (AoU)"
full_results$group4 <- "RICE-RV; Train (AoU) -> Validate (AoU)"
full_results$p.signif_beta2 <- ""
full_results$p.signif_beta2[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] > 0,"**",""))
full_results$position2 <- NA
full_results$position2[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] <- full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"] + 0.03


ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV; Train (AoU) -> Validate (UKB)"] + full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (UKB)"],
              full_results$beta_adjusted[full_results$Method == "RICE-CV; Train (AoU) -> Validate (AoU)"] + full_results$beta_adjusted[full_results$Method == "RICE-RV; Train (AoU) -> Validate (AoU)"],
              full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("RICE-CV; Train (AoU) -> Validate (AoU)","RICE-CV; Train (AoU) -> Validate (UKB)","RICE-RV; Train (AoU) -> Validate (AoU)","RICE-RV; Train (AoU) -> Validate (UKB)"))

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) +
  ylab("Beta of PRS per SD") +
  ylim(0,ylim) +
  stat_pvalue_manual(full_results,
                     label = "p.signif_beta1",
                     y.position = "position1",
                     xmin = "group1",
                     xmax = "group2",
                     size = 2.5) +
  stat_pvalue_manual(full_results,
                     label = "p.signif_beta2",
                     y.position = "position2",
                     xmin = "group3",
                     xmax = "group4",
                     size = 2.5) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="Method + Dataset"))

ggsave(paste0("AoU_CrossPlatform_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV; Train (AoU) -> Validate (AoU)",]
full_results <- full_results[full_results$Method != "RICE-RV; Train (AoU) -> Validate (UKB)",]
full_results$Method[full_results$Method == "RICE-CV; Train (AoU) -> Validate (AoU)"] <- "RICE; Train (AoU) -> Validate (AoU)"
full_results$Method[full_results$Method == "RICE-CV; Train (AoU) -> Validate (UKB)"] <- "RICE; Train (AoU) -> Validate (UKB)"
full_results$Method <- factor(full_results$Method,levels = c("RICE; Train (AoU) -> Validate (AoU)","RICE; Train (AoU) -> Validate (UKB)"))

ylim <- max(full_results$R2_adjusted) + 0.03

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method, y=abs(R2_adjusted),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) +
  ylab("R2") +
  ylim(0,ylim) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="Method + Dataset"))

ggsave(paste0("AoU_CrossPlatform_Adjusted_R2.png"),g2,width=10, height=6.18047,dpi = 300)