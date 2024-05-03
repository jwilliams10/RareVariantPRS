rm(list = ls())

WES_Results_Continuous <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Continuous.csv")
WES_Results_Continuous$Data_Type <- "WES"
WES_Results_Binary <- read.csv("~/Desktop/RareVariantPRS_Results/WES_Results_Binary.csv")
WES_Results_Binary$Data_Type <- "WES"
WGS_Results_Binary <- read.csv("~/Desktop/RareVariantPRS_Results/WGS_Results_Binary.csv")
WGS_Results_Binary$Data_Type <- "WGS"
WGS_Results_Continuous <- read.csv("~/Desktop/RareVariantPRS_Results/WGS_Results_Continuous.csv")
WGS_Results_Continuous$Data_Type <- "WGS"

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

full_results_Continuous <- rbind(WES_Results_Continuous,WGS_Results_Continuous)
full_results_Continuous <- full_results_Continuous[full_results_Continuous$Method %in% c("CV","RV"),]
full_results_Continuous$Method_DataSource <- paste0(full_results_Continuous$Method,"_",full_results_Continuous$Data_Type)

ggplot(full_results_Continuous) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(vars(trait), vars(ancestry)) +
  ggtitle("WES vs WGS Adjusted PRS Results") +
  ylab("Beta") +
  ylim(0,0.6) +
  theme_Publication() +
  scale_fill_Publication()