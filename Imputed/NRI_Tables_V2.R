rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)

continuous_traits <- c("Height","BMI","TC","HDL","LDL","logTG")
NRI_Data_Continuous <- NULL
for(trait in continuous_traits){
  
  pheno_validation <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  CV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble/",trait,"_PRS_Validation.csv"))
  colnames(CV_PRS_Validation) <- c("IID","CV_PRS")
  pheno_validation <- inner_join(pheno_validation,CV_PRS_Validation)
  RV_PRS_Validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/SingleTrait_Ensemble_RV/",trait,"_PRS_Validation.csv"))
  colnames(RV_PRS_Validation) <- c("IID","RV_PRS")
  pheno_validation <- inner_join(pheno_validation,RV_PRS_Validation)
  
  model.null <- lm(as.formula(paste0(trait,"~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data=pheno_validation)
  pheno_validation$y_validation <- NA
  pheno_validation$y_validation[!is.na(pheno_validation[,trait])] <- model.null$residual
  
  pheno_validation <- pheno_validation[!is.na(pheno_validation[,trait]),]
  
  CV_RV_PRS_raw <- pheno_validation
  CV_RV_PRS_adjusted <- pheno_validation
  
  for(i in c("RV_PRS","CV_PRS")){
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
  
  NRI_Data_Continuous <- NULL
  
  for(risk in c(0.05,0.1)){
    
    truth_HighRisk <- which(pheno_validation$y_validation > quantile(pheno_validation$y_validation,0.9))
    CV_HighRisk_RV_HighRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_HighRisk_RV_NotHighRisk <- which((pheno_validation$CV_PRS > quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,1 - risk)))
    
    truth_NotHighRisk <- which(pheno_validation$y_validation < quantile(pheno_validation$y_validation,0.9))
    CV_NotHighRisk_RV_HighRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS > quantile(pheno_validation$RV_PRS,1 - risk)))
    CV_NotHighRisk_RV_NotHighRisk <- which((pheno_validation$CV_PRS < quantile(pheno_validation$CV_PRS,0.9)) & (pheno_validation$RV_PRS < quantile(pheno_validation$RV_PRS,1 - risk)))
    
    tmp <- data.frame(trait = trait, risk = paste0(100*risk,"%"),
                      
                      Percent_A = 100*length(CV_NotHighRisk_RV_NotHighRisk)/nrow(pheno_validation),
                      Mean_A = mean(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk]),
                      SE_A = sd(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk])/sqrt(length(CV_NotHighRisk_RV_NotHighRisk)),
                      Total_Capture_A = sum(CV_NotHighRisk_RV_NotHighRisk %in% truth_HighRisk),
                      
                      Percent_B = 100*length(CV_HighRisk_RV_NotHighRisk)/nrow(pheno_validation),
                      Mean_B = mean(pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk]),
                      SE_B = sd(pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk])/sqrt(length(CV_HighRisk_RV_NotHighRisk)),
                      Total_Capture_B = sum(CV_HighRisk_RV_NotHighRisk %in% truth_HighRisk),
                      
                      Percent_C = 100*length(CV_NotHighRisk_RV_HighRisk)/nrow(pheno_validation),
                      Mean_C = mean(pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk]),
                      SE_C = sd(pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk])/sqrt(length(CV_NotHighRisk_RV_HighRisk)),
                      Total_Capture_C = sum(CV_NotHighRisk_RV_HighRisk %in% truth_HighRisk),
                      
                      Percent_D = 100*length(CV_HighRisk_RV_HighRisk)/nrow(pheno_validation),
                      Mean_D = mean(pheno_validation$y_validation[CV_HighRisk_RV_HighRisk]),
                      SE_D = sd(pheno_validation$y_validation[CV_HighRisk_RV_HighRisk])/sqrt(length(CV_HighRisk_RV_HighRisk)),
                      Total_Capture_D = sum(CV_HighRisk_RV_HighRisk %in% truth_HighRisk),
                      
                      Mean_BD = mean(pheno_validation$y_validation[c(CV_HighRisk_RV_HighRisk,CV_HighRisk_RV_NotHighRisk)]))
    
    tmp$C_minus_A <- (tmp$Mean_C - tmp$Mean_A)/sd(pheno_validation$y_validation)
    tmp$B_minus_A <- (tmp$Mean_B - tmp$Mean_A)/sd(pheno_validation$y_validation)
    tmp$BD_minus_A <- (tmp$Mean_BD - tmp$Mean_A)/sd(pheno_validation$y_validation)
    
    tmp$A_vs_B <- t.test(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk], pheno_validation$y_validation[CV_HighRisk_RV_NotHighRisk], alternative = "two.sided", var.equal = FALSE)$p.value
    tmp$A_vs_C <- t.test(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk], pheno_validation$y_validation[CV_NotHighRisk_RV_HighRisk], alternative = "two.sided", var.equal = FALSE)$p.value
    tmp$A_vs_D <- t.test(pheno_validation$y_validation[CV_NotHighRisk_RV_NotHighRisk], pheno_validation$y_validation[CV_HighRisk_RV_HighRisk], alternative = "two.sided", var.equal = FALSE)$p.value
    
    NRI_Data_Continuous <- rbind(NRI_Data_Continuous,tmp)
  } 
  
  CV_RV_PRS_adjusted$Common_Bin <- 9
  CV_RV_PRS_adjusted$Rare_Bin <- 9
  
  Common_quants <- quantile(CV_RV_PRS_adjusted$CV_PRS,c(0,.1,.2,.3,.4,.6,.7,.8,.9))
  
  for(i in 1:8){
    CV_RV_PRS_adjusted$Common_Bin[CV_RV_PRS_adjusted$CV_PRS < unname(Common_quants[i + 1]) & CV_RV_PRS_adjusted$CV_PRS >= unname(Common_quants[i])] <- i
  }
  
  Rare_quants <- quantile(CV_RV_PRS_adjusted$RV_PRS,c(0,.05,.2,.3,.4,.6,.7,.8,.95))
  for(i in 1:8){
    CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$RV_PRS < unname(Rare_quants[i + 1]) & CV_RV_PRS_adjusted$RV_PRS >= unname(Rare_quants[i])] <- i
  }
  
  CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$Rare_Bin == 1] <- "Below 5%"
  CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$Rare_Bin %in% c("4","5","6")] <- "30% - 70%"
  CV_RV_PRS_adjusted$Rare_Bin[CV_RV_PRS_adjusted$Rare_Bin %in% c("9")] <- "Above 95%"
  
  CV_RV_PRS_adjusted$y_validation <- scale(CV_RV_PRS_adjusted$y_validation)
  
  CV_RV_PRS_adjusted <- CV_RV_PRS_adjusted[CV_RV_PRS_adjusted$Rare_Bin %in% c("Below 5%","30% - 70%","Above 95%"),]
  
  CV_RV_PRS_adjusted$Rare_Bin <- factor(CV_RV_PRS_adjusted$Rare_Bin,levels = c("Below 5%","30% - 70%","Above 95%"))
  
  CV_RV_PRS_adjusted_se <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted,function(x){sd(x)/sqrt(length(x))})
  CV_RV_PRS_adjusted <- aggregate(y_validation ~ Common_Bin + Rare_Bin,data = CV_RV_PRS_adjusted,mean)
  
  colnames(CV_RV_PRS_adjusted) <- c("Common_Bin","Rare_Bin","Mean")
  colnames(CV_RV_PRS_adjusted_se) <- c("Common_Bin","Rare_Bin","SE")
  CV_RV_PRS_adjusted <- inner_join(CV_RV_PRS_adjusted,CV_RV_PRS_adjusted_se)
  
  colnames(CV_RV_PRS_adjusted) <- c("Common_Bin","RICE-RV Quantiles (Rare Variants)","Mean","SE")
  
  ymin <- round(min(c(CV_RV_PRS_adjusted$Mean - CV_RV_PRS_adjusted$SE)) - 0.05,2)
  ymax <- round(max(c(CV_RV_PRS_adjusted$Mean + CV_RV_PRS_adjusted$SE)) + 0.05,2)
  
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
              strip.text.x = element_text(face = "bold",size = 16),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.line = element_line(colour="black",size=2),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_line(),
              # panel.grid.major = element_line(colour="#f0f0f0"),
              # panel.grid.minor = element_line(colour="#f0f0f0"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.text=element_text(size=12),
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
  
  plot1 <- ggplot(data=CV_RV_PRS_adjusted, aes(x=Common_Bin, y=Mean, color=`RICE-RV Quantiles (Rare Variants)`)) + geom_line() + geom_pointrange(aes(ymin=Mean-SE, ymax=Mean+SE)) + theme_Publication() + ylab(paste0("Standardized ",ifelse(trait == "logTG","log(TG)",trait))) + ylim(c(ymin,ymax)) +
    scale_x_continuous(breaks = c(1:9),labels = c("0-10%","10-20%","20-30%","30-40%","40-60%","60-70%","70-80%","80-90%","90-100%")) + labs(x = "RICE-CV Quantiles (Common Variants)") + ggtitle(paste0("Standardized ",trait," Across RICE-CV Quantiles, Stratified by RICE-RV Quantiles"))
  
  
  
  
  plot_data <- data.frame(Method = rep(c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
                          trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
                          risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
                          value = c(NRI_Data_Continuous$Total_Capture_A,NRI_Data_Continuous$Total_Capture_B,NRI_Data_Continuous$Total_Capture_C,NRI_Data_Continuous$Total_Capture_D))
  
  plot_data$Method <- factor(plot_data$Method, levels = c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"))
  
  theme_Publication <- function(base_size=12) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, )
      + theme(plot.title = element_text(face = "bold",
                                        size = 16, hjust = 0.5),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              # panel.grid.major = element_line(colour="#f0f0f0"),
              # panel.grid.minor = element_line(colour="#f0f0f0"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              #legend.position = "bottom",
              #legend.direction = "horizontal",
              #legend.key.size= unit(0.2, "cm"),
              #legend.margin = unit(0, "cm"),
              legend.title = element_text(face="bold.italic", size = 16),
              legend.text = element_text(size = 12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold"),
              legend.position = "bottom",             # Place legend at bottom
              legend.direction = "horizontal" 
      ))
    
  }
  
  scale_fill_Publication <- function(...){
    library(scales)
    discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#973999","#009cdf")), ...)
  }
  
  scale_color_Publication <- function(...){
    library(scales)
    discrete_scale("color","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#973999","#009cdf")), ...)
  }
  
  plot_data_sub <- plot_data[plot_data$trait == trait & plot_data$risk == "5%", ] %>%
    arrange(Method) %>%
    mutate(
      fraction = value / sum(value),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, -1)),
      label_pos = (ymin + ymax) / 2,
      label = paste0(round(100*value/sum(value), 1), "%")
    )
  
  plot_data_sub$Method <- factor(plot_data_sub$Method, levels = c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"))
  
  plot2 <- ggplot(plot_data_sub, aes(ymax = ymax, ymin = ymin, xmax = 1.5, xmin = 0, fill = Method)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y") +
    geom_text_repel(
      aes(x = 1.5, y = label_pos, label = label),
      nudge_x = 0.3,     # push labels outward
      size = 5,
      fontface = "bold",
      segment.color = "black", # connector lines
      show.legend = FALSE
    ) + 
    guides(color = guide_legend(ncol = 2)) + 
    ggtitle(paste0("Proportion of High ",trait," Phenotype Individuals\n(Top 10% Quantile) Captured by RICE-CV\nand RICE-RV Risk Groups")) + 
    theme_Publication() +
    scale_fill_Publication()
  
  theme_Publication <- function(base_size=12) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, )
      + theme(plot.title = element_text(face = "bold",
                                        size = 16, hjust = 0.5),
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
              legend.title = element_text(face="bold.italic", size = 16),
              legend.text = element_text(size = 12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold"),
              legend.position = "bottom",             # Place legend at bottom
              legend.direction = "horizontal" 
      ))
    
  }
  
  plot_data <- data.frame(Method = rep(c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
                          trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
                          risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
                          mean = c(NRI_Data_Continuous$Mean_A,NRI_Data_Continuous$Mean_B,NRI_Data_Continuous$Mean_C,NRI_Data_Continuous$Mean_D),
                          se = c(NRI_Data_Continuous$SE_A,NRI_Data_Continuous$SE_B,NRI_Data_Continuous$SE_C,NRI_Data_Continuous$SE_D))
  
  plot_data$Method <- factor(plot_data$Method, levels = c("Low CV PRS, Low RV PRS","High CV PRS, Low RV PRS","Low CV PRS, High RV PRS","High CV PRS, High RV PRS"))
  
  if(trait %in% c("LDL","BMI","TC")){
    stat.test <- data.frame(
      group1 = "Low CV PRS, Low RV PRS",
      group2 = rep(c("High CV PRS, Low RV PRS", "Low CV PRS, High RV PRS", "High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
      trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
      risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
      p.adj = signif(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D),3),      # Adjusted p-values
      y.position = as.vector(outer(apply(cbind(NRI_Data_Continuous$Mean_B,NRI_Data_Continuous$Mean_C,NRI_Data_Continuous$Mean_D),1,max),c(1.3,1.4,1.5),"*")),     # Height of brackets
      p.adj.signif = ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) > 0.05, "",ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) < 0.01,"**","*"))  # Significance symbols
    ) 
  }else{
    stat.test <- data.frame(
      group1 = "Low CV PRS, Low RV PRS",
      group2 = rep(c("High CV PRS, Low RV PRS", "Low CV PRS, High RV PRS", "High CV PRS, High RV PRS"),each = nrow(NRI_Data_Continuous)),
      trait = c(NRI_Data_Continuous$trait,NRI_Data_Continuous$trait,NRI_Data_Continuous$trait),
      risk = c(NRI_Data_Continuous$risk,NRI_Data_Continuous$risk,NRI_Data_Continuous$risk),
      p.adj = signif(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D),3),      # Adjusted p-values
      y.position = as.vector(outer(apply(cbind(NRI_Data_Continuous$Mean_B,NRI_Data_Continuous$Mean_C,NRI_Data_Continuous$Mean_D),1,max),c(1.1,1.2,1.3),"*")),     # Height of brackets
      p.adj.signif = ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) > 0.05, "",ifelse(c(NRI_Data_Continuous$A_vs_B, NRI_Data_Continuous$A_vs_C, NRI_Data_Continuous$A_vs_D) < 0.01,"**","*"))  # Significance symbols
    )
  }
  
  # Plot
  plot3 <- ggplot(plot_data[plot_data$trait == trait & plot_data$risk == "5%",], aes(x = Method, y = mean,color = Method)) +
    geom_point(size = 4.5) +
    geom_linerange(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se)) +
    stat_pvalue_manual(stat.test[stat.test$trait == trait & stat.test$risk == "5%",], label = "p.adj", tip.length = 0.01) +
    theme_Publication() + 
    scale_color_Publication() +
    ggtitle(paste0("Mean Standardized ",trait," Levels by RICE-CV\nand RICE-RV Risk Groups")) + 
    labs(y = "Mean (95% CI)", x = "Group") 
  
  
  legend1 <- ggplotGrob(plot1)$grobs[[which(sapply(ggplotGrob(plot1)$grobs, function(x) x$name) == "guide-box")]]
  row1 <- plot_grid(plot1 + theme(legend.position="none"), legend1, ncol = 1, rel_heights = c(1, .15),labels = c("a",NULL),label_size = 18)
  legend2 <- ggplotGrob(plot2)$grobs[[which(sapply(ggplotGrob(plot2)$grobs, function(x) x$name) == "guide-box")]]
  row2 <- plot_grid(plot2 + theme(legend.position="none"), plot3 + theme(legend.position="none"), ncol = 2, rel_widths = c(1, 1),labels = c("b","c"),label_size = 18)
  row2 <- plot_grid(row2, legend2, ncol = 1, rel_heights = c(1, .15))
  
  final_plot <- plot_grid(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"),plot3 + theme(legend.position="none"), nrow = 1, rel_widths = c(1, 1, 1))
  
  pdf(paste0(trait,"_Risk.pdf"), width=12, height=12)
  
  print(plot_grid(row1,row2, ncol = 1, rel_heights = c(1.1, .9)))
  
  dev.off()
}