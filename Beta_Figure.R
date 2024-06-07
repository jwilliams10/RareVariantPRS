# rm(list = ls())
# 
# library(ggplot2)
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
#             panel.border = element_rect(colour = NA),
#             axis.title = element_text(face = "bold",size = 24),
#             axis.title.y = element_text(angle=90,vjust =2),
#             axis.text.x = element_blank(),
#             axis.text.y = element_text(size = 18),
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
# scale_fill_Publication <- function(...){
#   library(scales)
#   discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#E23838","#973999","#009cdf")), ...)
#   
# }
# 
# scale_color_Publication <- function(...){
#   library(scales)
#   discrete_scale("color","Publication",manual_pal(values = c("#973999","#E23838","#F78200","#E23838","#973999","#009cdf")), ...)
#   
# }
# 
# 
# data <- data.frame(height = rep(1,6),x_value = c("1","2","2","3","3","3"),color = c("1"),labels = c("Beta of PRS per SD","Beta of PRS per SD","Beta of PRS per SD"))
# data$x_value <- factor(data$x_value,levels = c("3","2","1"))
# 
# ggplot(data,aes(x=x_value, y=height,label = labels)) + geom_bar(stat="identity", alpha=0.7,fill="white",color = "#973999",size=2) + 
#   theme_Publication() + 
#   scale_fill_Publication() + theme(legend.position="none") + ylab("Expected Increase in Standardized Trait") + xlab("Change in SD of Standardized RICE-CV PRS") +
#   geom_text(size = 6, position = position_stack(vjust = 0.5),color = "black") + 
#   coord_flip() 
# 
# 
# data <- data.frame(height = rep(1,6),x_value = c("1","2","2","3","3","3"),color = c("1"),labels = c("Beta of PRS per SD","Beta of PRS per SD","Beta of PRS per SD"))
# data$x_value <- factor(data$x_value,levels = c("3","2","1"))
# 
# ggplot(data,aes(x=x_value, y=height,label = labels)) + geom_bar(stat="identity", alpha=0.7,fill="white",color = "#E23838",size=2) + 
#   theme_Publication() + 
#   scale_fill_Publication() + theme(legend.position="none") + ylab("Expected Increase in Standardized Trait") + xlab("Change in SD of Standardized RICE-RV PRS") +
#   geom_text(size = 6, position = position_stack(vjust = 0.5),color = "black") + 
#   coord_flip() 
# 
# scale_fill_Publication <- function(...){
#   library(scales)
#   discrete_scale("fill","Publication",manual_pal(values = c("#E23838","#FFB900","#F78200","#E23838","#973999","#009cdf")), ...)
#   
# }
# 
# ggplot(data,aes(x=x_value, y=height,fill = color,label = labels)) + geom_bar(stat="identity", alpha=0.7,color = "black",size=2) + 
#   theme_Publication() + 
#   scale_fill_Publication() + theme(legend.position="none") + ylab("Expected Increase in Standardized Trait") + xlab("Change in SD of Standardized RICE-RV PRS") +
#   geom_text(size = 6, position = position_stack(vjust = 0.5),color = "black") + 
#   coord_flip() 
# 
# ggplot(data,aes(x=x_value, y=height,label = labels)) + geom_bar(stat="identity", alpha=0.7,color = "black",size=2,fill = "white") + 
#   theme_Publication() + 
#   scale_fill_Publication() + theme(legend.position="none") + ylab("Expected Increase in Standardized Trait") + xlab("Change in SD of Standardized RICE-RV PRS") +
#   geom_text(size = 6, position = position_stack(vjust = 0.5),aes(color = color)) + 
#   coord_flip() 
# 
# 
# data <- data.frame(height = rep(1,12),x_value = c("1","1","2","2","2","2","3","3","3","3","3","3"),color = c("1","2","1","1","2","2","1","1","1","2","2","2"),labels = c("Beta of PRS/SD","Beta of PRS/SD","Beta of PRS/SD","Beta of PRS/SD","Beta of PRS/SD","Beta of PRS/SD"))
# data$x_value <- factor(data$x_value,levels = c("3","2","1"))
# 
# ggplot(data,aes(x=x_value, y=height,colour = color,label = labels)) + geom_bar(stat="identity", alpha=0.7,fill="white",size=2) + 
#   theme_Publication() + 
#   scale_color_Publication() + theme(legend.position="none") + ylab("Expected Increase in Standardized Trait") + xlab("Change in SD of Standardized RICE-RV PRS") +
#   geom_text(size = 6, position = position_stack(vjust = 0.5),color = "black") + 
#   coord_flip() 














































rm(list = ls())

library(ggplot2)
library(cowplot)

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
            axis.title = element_text(face = "bold",size = 24),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(face = "bold",size = 18),
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
            legend.title = element_blank(),
            legend.text = element_text(face ="bold",size = 18),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#E23838","#973999","#F78200","#E23838","#973999","#009cdf")), ...)
  
}

data <- data.frame(height = c(1,.5,1,.5),x_value = c("2","1","3","3"),Method = c("RICE-CV","RICE-RV","RICE-CV","RICE-RV"))
data$Method <- factor(data$Method,levels = c("RICE-RV","RICE-CV"))

ggplot(data) + geom_bar(aes(x = x_value, y = height, fill = Method), stat = "identity", alpha=0.7,width = 0.7) + 
  theme_Publication() + 
  scale_fill_Publication(labels = c("+1 SD of RICE-RV PRS","+1 SD of RICE-CV PRS")) +
  ylab("Change in Standardized Trait") + 
  scale_y_continuous(breaks = c(.5,1,1.5),labels = c(expression(beta[RV]),expression(beta[CV]),expression(beta[CV + RV])))





