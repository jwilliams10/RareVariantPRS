rm(list = ls())

library(ggplot2)
library(dplyr)

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            # panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
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

h2_dat <- read.csv("/data/williamsjacr/UKB_WES_Simulation/Sim_Characteristics.csv")

h2_dat_compressed <- aggregate(.~Ancestry + causal_prop + scaled,data = h2_dat[,c("Ancestry","causal_prop","scaled","h2_rare","Average_Burden_rare")],mean)
colnames(h2_dat_compressed) <- c("Ancestry","Causal_Prop","Scaled","Average_h2","Average_Burden")
h2_dat_compressed_se <- aggregate(h2_rare~Ancestry + causal_prop + scaled,data = h2_dat,function(x){quantile(x,0.025)})
colnames(h2_dat_compressed_se) <- c("Ancestry","Causal_Prop","Scaled","Q_025")
h2_dat_compressed <- inner_join(h2_dat_compressed,h2_dat_compressed_se)
h2_dat_compressed_se <- aggregate(h2_rare~Ancestry + causal_prop + scaled,data = h2_dat,function(x){quantile(x,0.975)})
colnames(h2_dat_compressed_se) <- c("Ancestry","Causal_Prop","Scaled","Q_975")
h2_dat_compressed <- inner_join(h2_dat_compressed,h2_dat_compressed_se)


h2_dat_compressed$Causal_Prop <- paste0("Causal Prop. ",h2_dat_compressed$Causal_Prop)
h2_dat_compressed$Scaled <- paste0("Scaled: ",h2_dat_compressed$Scaled)
h2_dat_compressed <- h2_dat_compressed[h2_dat_compressed$Ancestry != "EAS",]
h2_dat_compressed <- h2_dat_compressed[h2_dat_compressed$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01"),]
# Just first plot, similar as quantile plots, points +- se

g1 <- ggplot(data=h2_dat_compressed, aes(x=Average_Burden, y=Average_h2, color=Ancestry)) + geom_pointrange(aes(ymin=Q_025, ymax=Q_975)) + 
  theme_Publication() + ylab("Average Heritability") + labs(x = "Average of Causal Burden Scores") + facet_grid(Causal_Prop~Scaled)

ggsave(paste0("Sim_Characteristics.png"),g1,width=10, height=6.18047,dpi = 300)

# h2_dat$causal_prop <- paste0("Causal Prop. ",h2_dat$causal_prop)
# h2_dat$scaled <- paste0("Scaled: ",h2_dat$scaled)
# h2_dat <- h2_dat[h2_dat$Ancestry != "EAS",]
# h2_dat <- h2_dat[h2_dat$causal_prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01"),]
# 
# 
# theme_Publication <- function(base_size=12) {
#   library(grid)
#   library(ggthemes)
#   (theme_foundation(base_size=base_size, )
#     + theme(plot.title = element_text(face = "bold",
#                                       size = rel(1.1), hjust = 0.5),
#             text = element_text(),
#             panel.background = element_rect(colour = NA),
#             plot.background = element_rect(colour = NA),
#             # panel.border = element_rect(colour = NA),
#             axis.title = element_text(face = "bold",size = 16),
#             axis.title.y = element_text(angle=90,vjust =2),
#             axis.text.x = element_text(size = 10), 
#             axis.text.y = element_text(size = 10), 
#             axis.line = element_line(colour="black",size=2),
#             axis.ticks = element_line(),
#             # panel.grid.major = element_line(colour="#f0f0f0"),
#             # panel.grid.minor = element_line(colour="#f0f0f0"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             legend.key = element_rect(colour = NA),
#             #legend.position = "bottom",
#             #legend.direction = "horizontal",
#             #legend.key.size= unit(0.2, "cm"),
#             #legend.margin = unit(0, "cm"),
#             legend.title = element_text(face="bold.italic", size =18),
#             #legend.text = element_text(face ="bold"),
#             plot.margin=unit(c(10,5,5,5),"mm"),
#             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#             strip.text = element_text(face="bold")
#     ))
#   
# }
# 
# ggplot(data=h2_dat, aes(x=Average_Burden_rare, y=h2_rare, color=Ancestry,group = Average_Burden_rare)) + geom_boxplot() + 
#   theme_Publication() + ylab("Average Heritability") + labs(x = "Average of Causal Burden Scores") + facet_grid(causal_prop~scaled)


# ggplot(h2_dat) + geom_point(aes(x = Var_Burden_rare,y = h2_rare,color = Ancestry)) + facet_grid(causal_prop~scaled)
# 
# ggplot(h2_dat) + geom_point(aes(x = MAF_common,y = h2_common,color = Ancestry)) + facet_grid(causal_prop~scaled)
# ggplot(h2_dat) + geom_point(aes(x = Var_common,y = h2_common,color = Ancestry)) + facet_grid(causal_prop~scaled)
