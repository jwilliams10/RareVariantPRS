rm(list = ls())

library(ggplot2)

trait <- "LDL"

RV_PRS <- read.csv(paste0("~/Desktop/RareVariantPRS_Results/PRS_Validation/WGS/",trait,"_BestPRS.csv"))
CV_PRS <- read.delim(paste0("~/Desktop/RareVariantPRS_Results/PRS_Validation/WGS/",trait,"_Best_Validation_All.txt"))

load("~/Desktop/RareVariantPRS_Results/all_phenotypes.RData")

CV_RV_PRS <- inner_join(RV_PRS,CV_PRS)
CV_RV_PRS_adjusted <- CV_RV_PRS

for(i in c("RV_PRS","prs")){
  tmp <- data.frame(y = CV_RV_PRS_adjusted[,i],CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
  mod <- lm(y~.,data = tmp)
  R <- mod$residuals
  tmp <- data.frame(y = R^2,CV_RV_PRS_adjusted[,c("pc1","pc2","pc3","pc4","pc5")])
  mod <- lm(y~.,data = tmp)
  y_hat <- predict(mod,tmp)
  if(sum(y_hat < 0) > 0){
    mod <- lm(y~1,data = tmp)
    y_hat <- predict(mod,tmp)
  }
  if(sum(sqrt(y_hat)) == 0){
    CV_RV_PRS_adjusted[,i] <- 0
  }else{
    CV_RV_PRS_adjusted[,i] <- R/sqrt(y_hat)
  }
}

CV_RV_PRS_adjusted$prs <- scale(CV_RV_PRS_adjusted$prs)
CV_RV_PRS_adjusted$RV_PRS <- scale(CV_RV_PRS_adjusted$RV_PRS)

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

risk_cv <- quantile(CV_RV_PRS_adjusted$prs,c(.90))
risk_rv <- quantile(CV_RV_PRS_adjusted$RV_PRS[abs(CV_RV_PRS_adjusted$RV_PRS) > 0.01 & abs(CV_RV_PRS_adjusted$RV_PRS) < 7],c(.9))

p1 <- ggplot(CV_RV_PRS_adjusted, aes(x=prs)) + 
  geom_histogram(color="black", fill="#973999",aes(y=..density..),bins = 100) + 
  xlab("Standardized RICE-CV PRS") + ylab(element_blank()) + 
  theme_Publication() + scale_x_continuous(breaks = c(-5:5)) + 
  geom_vline(xintercept = as.numeric(risk_cv),color = "black",linetype = "dashed",size = 1.2) + 
  annotate("text", x=2.2, y=0.4, label= "Upper 10% Quantile",size = 5) + theme(text = element_text(size=20))

p1 
ggsave(p1,filename="Desktop/RareVariantPRS_Results/Figures/Fig2_CV.pdf",width = 10,height = 6.1804697157)

p2 <- ggplot(CV_RV_PRS_adjusted[abs(CV_RV_PRS_adjusted$RV_PRS) > 0.01 & abs(CV_RV_PRS_adjusted$RV_PRS) < 7,], aes(x=RV_PRS)) + 
  geom_histogram(color="black", fill="#E23838",aes(y=..density..),bins = 100) + 
  xlab("Standardized RICE-RV PRS") + ylab(element_blank()) + 
  theme_Publication() + scale_x_continuous(breaks = c(-7:7)) + 
  geom_vline(xintercept = as.numeric(risk_rv),color = "black",linetype = "dashed",size = 1.2) + 
  annotate("text", x=2.4, y=2.5, label= "Upper 10% Quantile",size = 5) + theme(text = element_text(size=20))

p2
ggsave(p2,filename="Desktop/RareVariantPRS_Results/Figures/Fig2_RV.pdf",width = 10,height = 6.1804697157)

# tmp <- inner_join(ukb_pheno[,c("IID",trait)],CV_RV_PRS_adjusted)
# plot(tmp[,"RV_PRS"],tmp[,trait])



