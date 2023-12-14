rm(list = ls())

library(ggplot2)

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
            axis.title.x = element_text(vjust = -0.2),
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

a <- data.frame(Trait = NULL,Ancestry = NULL,Method = NULL,AUC = NULL,AUC_Low = NULL,AUC_High = NULL)

Trait <- "Asthma"
Ancestry <- "NonEur"

for(Trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  for(Ancestry in c("Eur","NonEur","UNK","MIX","AFR","SAS","EAS")){
    
    if(Ancestry == "Eur"){
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",Trait,"_CT_result_",toupper(Ancestry),".RData")) 
    }else{
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",Trait,"_CT_result_",Ancestry,".RData"))
    }
    if(sum(is.na(ct.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "CT",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "CT",AUC = ct.result[[1]]$AUC,AUC_Low = ct.result[[1]]$AUC_low,AUC_High = ct.result[[1]]$AUC_high))
    }
    
    if(Ancestry == "Eur"){
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",Trait,"_ldpred2_result_",toupper(Ancestry),".RData"))
    }else{
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",Trait,"_ldpred2_result_",Ancestry,".RData"))
    }
    if(sum(is.na(ldpred2.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LDPred2",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LDPred2",AUC = ldpred2.result[[1]]$AUC,AUC_Low = ldpred2.result[[1]]$AUC_low,AUC_High = ldpred2.result[[1]]$AUC_high))
    }
    
    if(Ancestry == "Eur"){
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",Trait,"_LASSOSUM2_result_",toupper(Ancestry),".RData"))
    }else{
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",Trait,"_LASSOSUM2_result_",Ancestry,".RData"))
    }
    if(sum(is.na(ldpred2.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LASSOSUM2",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LASSOSUM2",AUC = LASSOSUM2.result[[1]]$AUC,AUC_Low = LASSOSUM2.result[[1]]$AUC_low,AUC_High = LASSOSUM2.result[[1]]$AUC_high))
    }
    
    
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",Trait,"_sl_result_All_",Ancestry,".RData"))
    if(sum(is.na(SL.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "SL_Common",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "SL_Common",AUC = SL.result$AUC,AUC_Low = SL.result$AUC_low,AUC_High = SL.result$AUC_high))
    }
    
  }
}

results <- a
rm(a)

results <- results[!(results$Ancestry %in% c("EAS","MIX","UNK")),]

ggplot(results) +
  geom_bar(aes(x=Method, y=AUC,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=AUC_Low, ymax=AUC_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Common Variants") + 
  coord_cartesian(ylim=c(0.5,1)) + 
  ylab("AUC") + 
  theme_Publication() + 
  scale_fill_Publication()
