rm(list = ls())

library(ggplot2)

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

full_results_Continuous <- rbind(WES_Results_Continuous,WGS_Results_Continuous)
full_results_Continuous$Method[full_results_Continuous$Method == "CV"] <- "RICE-CV"
full_results_Continuous$Method[full_results_Continuous$Method == "RV"] <- "RICE-RV"
full_results_Continuous <- full_results_Continuous[full_results_Continuous$Method %in% c("RICE-CV","RICE-RV"),]
full_results_Continuous$Method_DataSource <- paste0(full_results_Continuous$Method," & ",full_results_Continuous$Data_Type)

full_results_Continuous <- full_results_Continuous[full_results_Continuous$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_Continuous$trait[full_results_Continuous$trait == "logTG"] <- "log(TG)"
full_results_Continuous$trait <- factor(full_results_Continuous$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))

full_results_Binary <- rbind(WES_Results_Binary,WGS_Results_Binary)
full_results_Binary$Method[full_results_Binary$Method == "CV"] <- "RICE-CV"
full_results_Binary$Method[full_results_Binary$Method == "RV"] <- "RICE-RV"
full_results_Binary <- full_results_Binary[full_results_Binary$Method %in% c("RICE-CV","RICE-RV"),]
full_results_Binary$Method_DataSource <- paste0(full_results_Binary$Method," & ",full_results_Binary$Data_Type)

full_results_Binary <- full_results_Binary[full_results_Binary$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_Binary$trait <- factor(full_results_Binary$trait,levels = c("Asthma","Breast","CAD","Prostate","T2D"))

plot1 <- ggplot(full_results_Continuous[full_results_Continuous$ancestry == "EUR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.7) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

plot2 <- ggplot(full_results_Binary[full_results_Binary$ancestry == "EUR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.95) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

prow <- plot_grid(NULL,
  plot1 + theme(legend.position="none") + ggtitle("Comparison of RICE-CV and RICE-RV PRS Results Using WES and WGS Data for European Ancestry"),
  NULL,
  plot2 + theme(legend.position="none") + theme(plot.title =element_blank()),
  rel_heights = c(-.05, 1, -0.05, 1),
  ncol = 1
)

legend_b <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]

pdf(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_vs_WGS_EUR.pdf"), width=10, height=6.18047)

print(plot_grid(prow,NULL, legend_b, ncol = 1, rel_heights = c(1,-.05, .1)))

dev.off()

plot1 <- ggplot(full_results_Continuous[full_results_Continuous$ancestry == "AFR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.7) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

plot2 <- ggplot(full_results_Binary[full_results_Binary$ancestry == "AFR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.95) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

prow <- plot_grid(NULL,
                  plot1 + theme(legend.position="none") + ggtitle("Comparison of RICE-CV and RICE-RV PRS Results Using WES and WGS Data for African Ancestry"),
                  NULL,
                  plot2 + theme(legend.position="none") + theme(plot.title =element_blank()),
                  rel_heights = c(-.05, 1, -0.05, 1),
                  ncol = 1
)

legend_b <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]

pdf(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_vs_WGS_AFR.pdf"), width=10, height=6.18047)

print(plot_grid(prow,NULL, legend_b, ncol = 1, rel_heights = c(1,-.05, .1)))

dev.off()

plot1 <- ggplot(full_results_Continuous[full_results_Continuous$ancestry == "AMR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.7) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

plot2 <- ggplot(full_results_Binary[full_results_Binary$ancestry == "AMR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.95) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

prow <- plot_grid(NULL,
                  plot1 + theme(legend.position="none") + ggtitle("Comparison of RICE-CV and RICE-RV PRS Results Using WES and WGS Data for Admixed American Ancestry"),
                  NULL,
                  plot2 + theme(legend.position="none") + theme(plot.title =element_blank()),
                  rel_heights = c(-.05, 1, -0.05, 1),
                  ncol = 1
)

legend_b <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]

pdf(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_vs_WGS_AMR.pdf"), width=10, height=6.18047)

print(plot_grid(prow,NULL, legend_b, ncol = 1, rel_heights = c(1,-.05, .1)))

dev.off()

plot1 <- ggplot(full_results_Continuous[full_results_Continuous$ancestry == "SAS",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.7) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

plot2 <- ggplot(full_results_Binary[full_results_Binary$ancestry == "SAS",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted),fill=Method_DataSource),position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0,0.95) +
  theme_Publication() +
  scale_fill_Publication() + guides(fill=guide_legend(title="PRS Method & Dataset"))

prow <- plot_grid(NULL,
                  plot1 + theme(legend.position="none") + ggtitle("Comparison of RICE-CV and RICE-RV PRS Results Using WES and WGS Data for South Asian Ancestry"),
                  NULL,
                  plot2 + theme(legend.position="none") + theme(plot.title =element_blank()),
                  rel_heights = c(-.05, 1, -0.05, 1),
                  ncol = 1
)

legend_b <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]

pdf(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_WES_vs_WGS_SAS.pdf"), width=10, height=6.18047)

print(plot_grid(prow,NULL, legend_b, ncol = 1, rel_heights = c(1,-.05, .1)))

dev.off()
