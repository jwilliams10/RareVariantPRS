rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(ggExtra)

trait <- "HDL"

pheno_validation <- read.delim("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/All_Validation.txt")
CV_PRS_Validation <- read.delim(paste0("/data/williamsjacr/UKB_WGS_Results/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestRareVariantPRS/",trait,"Validation_BestPRS.csv"))
colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)

model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
pheno_validation$y_validation <- NA
pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual

load("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/all_phenotypes.RData")

CV_RV_PRS <- pheno_validation
CV_RV_PRS_adjusted <- CV_RV_PRS

RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"Best_Betas.csv"))
RICE_CV_Results$Method <- "RICE-CV"
RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"Best_Betas.csv"))
RICE_RV_Results$Method <- "RICE-RV"

beta1 <- RICE_CV_Results$beta_adjusted[RICE_CV_Results$ancestry == "EUR" & RICE_CV_Results$trait == trait]

beta2 <- RICE_RV_Results$beta_adjusted[RICE_RV_Results$ancestry == "EUR" & RICE_RV_Results$trait == trait]

for(i in c("CV_PRS","RV_PRS")){
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
            axis.title.y = element_text(angle=90,vjust =2),
            axis.text.y = element_text(size = 20),
            axis.text.x = element_text(size = 20),
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

CV_RV_PRS_adjusted <- inner_join(CV_RV_PRS_adjusted[,c("IID","CV_PRS","RV_PRS")],ukb_pheno)
CV_RV_PRS_adjusted <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$ancestry == "EUR",]
CV_RV_PRS_adjusted[,trait] <- scale(CV_RV_PRS_adjusted[,trait])

risk_cv <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(.90))
risk_rv <- quantile(CV_RV_PRS_adjusted$RV_PRS,c(.90))

g <- ggplot(CV_RV_PRS_adjusted, aes_string(x="CV_PRS", y=trait)) + geom_point(alpha = .2) + 
  theme_Publication() + xlab("Standardized RICE-CV PRS") + 
  ylab(paste0("Standardized ",trait)) + geom_abline(intercept = 0, slope = beta1,col = "#973999",size = 1.5) +
  geom_vline(xintercept = as.numeric(risk_cv),color = "#306FBB",linetype = "dashed",size = 1.5) 
  # annotate("text", x=3, y=6, label= bquote(beta[1] == .(round(beta1,3))),size = 16/.pt)
p3 <- ggExtra::ggMarginal(g, type = "histogram",
                    xparams = list(color="black", fill="#973999",bins = 100),
                    yparams = list(color="black", fill="white",bins = 100))
p3

ggsave(p3,filename="Fig2_CV.pdf",width = 10,height = 10)

g <- ggplot(CV_RV_PRS_adjusted, aes_string(x="RV_PRS", y=trait)) + geom_point(alpha = .2) + 
  theme_Publication() + xlab("Standardized RICE-RV PRS") + 
  ylab(paste0("Standardized ",trait)) + geom_abline(intercept = 0, slope = beta2,col = "#E23838",size = 1.5) +
  geom_vline(xintercept = as.numeric(risk_rv),color = "#306FBB",linetype = "dashed",size = 1.5)
  # annotate("text", x=6, y=6, label= bquote(beta[2] == .(round(beta2,3))),size = 16/.pt)
p4 <- ggExtra::ggMarginal(g, type = "histogram",
                    xparams = list(color="black", fill="#E23838",bins = 100),
                    yparams = list(color="black", fill="white",bins = 100))
p4

ggsave(p4,filename="Fig2_RV.pdf",width = 10,height = 10)


p1 <- ggplot(CV_RV_PRS_adjusted, aes(x=CV_PRS)) +
  geom_histogram(color="black", fill="#973999",aes(y=..density..),bins = 100) +
  xlab("Standardized RICE-CV PRS") + ylab(element_blank()) +
  theme_Publication() + scale_x_continuous(breaks = c(-5:5)) +
  geom_vline(xintercept = as.numeric(risk_cv),color = "#5EBD3E",linetype = "dashed",size = 1.2)

p1
# ggsave(p1,filename="Desktop/RareVariantPRS_Results/Figures/Fig2_CV.pdf",width = 10,height = 6.1804697157)

p2 <- ggplot(CV_RV_PRS_adjusted[abs(CV_RV_PRS_adjusted$RV_PRS) > 0.05 & abs(CV_RV_PRS_adjusted$RV_PRS) < 7,], aes(x=RV_PRS)) +
  geom_histogram(color="black", fill="#E23838",aes(y=..density..),bins = 100) +
  xlab("Standardized RICE-RV PRS") + ylab(element_blank()) +
  theme_Publication() + scale_x_continuous(breaks = c(-7:7)) +
  geom_vline(xintercept = as.numeric(risk_rv),color = "#5EBD3E",linetype = "dashed",size = 1.2)

p2
# ggsave(p2,filename="Desktop/RareVariantPRS_Results/Figures/Fig2_RV.pdf",width = 10,height = 6.1804697157)

overlap_data <- CV_RV_PRS_adjusted
overlap_data$Final_PRS <- beta1*overlap_data$CV_PRS + beta2*overlap_data$RV_PRS
overlap_data$Variable <- ifelse(overlap_data$RV_PRS > 0.1,"0.1 < PRS",ifelse(overlap_data$RV_PRS < -0.1,"PRS < -0.1","-0.1 \u2264 PRS \u2264 0.1"))
overlap_data$Variable <- factor(overlap_data$Variable,levels = c("PRS < -0.1","-0.1 \u2264 PRS \u2264 0.1","0.1 < PRS"))
aggregate(Final_PRS~Variable,data = overlap_data,mean)

# overlap_data <- data.frame(Final_PRS = c(overlap_data$CV_PRS,overlap_data$Final_PRS,overlap_data[,trait]),Variable = rep(c("RICE-CV","RICE","Trait"),each = nrow(overlap_data)))

p3 <- ggplot(overlap_data,aes(x=Final_PRS, color=Variable)) + geom_density(alpha=0.25) + scale_color_manual(values = c("#973999","#5EBD3E","#E23838")) +
  theme_Publication() + guides(color=guide_legend(title="RICE-RV PRS")) + xlab("RICE PRS") + ylab(element_blank())
p3
# tmp <- inner_join(ukb_pheno[,c("IID",trait)],CV_RV_PRS_adjusted)
# plot(tmp[,"RV_PRS"],tmp[,trait])



