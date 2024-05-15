rm(list = ls())

# full_results <- NULL
# 
# for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
#   CT_Results <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"_Best_Betas_CT.csv"))
#   CT_Results$Method <- "CT"
#   LDPred2_Results <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"Best_Betas_LDPred2.csv"))
#   LDPred2_Results$Method <- "LDPred"
#   LASSOSUM2_Results <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"Best_Betas_LASSOSum.csv"))
#   LASSOSUM2_Results$Method <- "LASSOSum"
#   CV_Results <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"_Best_Betas_CV_SL.csv"))
#   CV_Results$Method <- "CV_SL"
#   
#   RV_Results_Coding <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"_Coding_Best_Betas.csv"))
#   RV_Results_Coding$Method <- "Coding"
#   
#   RV_Results_Noncoding <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"_Noncoding_Best_Betas.csv"))
#   RV_Results_Noncoding$Method <- "Noncoding"
#   
#   RV_Results_SL <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"_Best_Betas_RV_SL.csv"))
#   RV_Results_SL$Method <- "RV_SL"
#   
#   CV_RV_Results <- read.csv(paste0("/Users/williamsjacr/Downloads/Results/",trait,"_Best_Betas_CV_RV.csv"))
#   
#   
#   
#   full_results <- rbind(full_results,rbind(CT_Results,LDPred2_Results,LASSOSUM2_Results,CV_Results,RV_Results_Coding,RV_Results_Noncoding,RV_Results_SL,CV_RV_Results))
# }
# 
# full_results$beta_adjusted[full_results$beta_adjusted < 0 & full_results$Method %in% c("LDPred","LASSOSum")] <- -1*full_results$beta_adjusted[full_results$beta_adjusted < 0 & full_results$Method %in% c("LDPred","LASSOSum")]
# full_results$beta_raw[full_results$beta_raw < 0 & full_results$Method %in% c("LDPred","LASSOSum")] <- -1*full_results$beta_raw[full_results$beta_raw < 0 & full_results$Method %in% c("LDPred","LASSOSum")]
# 
# full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
# full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results <- read.csv("~/Desktop/RareVariantPRS_Results/WGS_Results_Continuous.csv")

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


library(ggplot2)
CV_Coding_Noncoding <- full_results[full_results$Method %in% c("Coding","Noncoding","CV","RV"),]
full_results <- full_results[!(full_results$Method %in% c("Coding","Noncoding","RV_SL")),]

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT","LDPred","LASSOSum","CV_SL","RV","CV"))
full_results$Method1[full_results$Method1 == "RV"] <- "CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT","LDPred","LASSOSum","CV_SL","CV"))

ggplot(full_results) +
  geom_bar(aes(x=Method1, y=abs(beta_raw),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("WGS Raw PRS Results") + 
  ylab("Beta") + 
  ylim(0,0.6) +
  theme_Publication() + 
  scale_fill_Publication()

ggplot(full_results) +
  geom_bar(aes(x=Method1, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("WGS Adjusted PRS Results") + 
  ylab("Beta") + 
  ylim(0,0.6) +
  theme_Publication() + 
  scale_fill_Publication()

ggplot(CV_Coding_Noncoding) +
  geom_bar(aes(x=Method, y=abs(beta_raw),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("WGS Raw PRS Results") + 
  ylab("Beta") + 
  ylim(0,0.6) +
  theme_Publication() + 
  scale_fill_Publication()
