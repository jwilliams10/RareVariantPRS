rm(list = ls())
library(ggplot2)
library(dplyr)

trait <- "BMI"

results <- NULL

for(trait in c("BMI","TC","HDL","LDL","logTG","Height")){
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EUR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EUR",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_NonEur_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "NonEur",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_AFR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "AFR",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_SAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "SAS",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EAS",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_UNK_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "UNK",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_MIX_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "MIX",ct.result[[1]]))
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EUR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EUR",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_NonEur.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "NonEur",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_AFR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "AFR",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_SAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "SAS",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EAS",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_UNK.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "UNK",ct.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_MIX.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "MIX",ct.result[[1]]))
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EUR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EUR",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_NonEur.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "NonEur",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_AFR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "AFR",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_SAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "SAS",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EAS",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_UNK.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "UNK",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_MIX.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "MIX",ldpred2.result[[1]]))
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EUR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EUR",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_NonEur_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "NonEur",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_AFR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "AFR",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_SAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "SAS",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EAS",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_UNK_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "UNK",ldpred2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_MIX_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "MIX",ldpred2.result[[1]]))
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EUR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EUR",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_NonEur.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "NonEur",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_AFR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "AFR",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_SAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "SAS",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EAS",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_UNK.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "UNK",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_MIX.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "MIX",LASSOSUM2.result[[1]]))
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EUR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EUR",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_NonEur_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "NonEur",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_AFR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "AFR",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_SAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "SAS",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EAS",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_UNK_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "UNK",LASSOSUM2.result[[1]]))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_MIX_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "MIX",LASSOSUM2.result[[1]]))
  
  
  
  
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EUR",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "NonEur",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "AFR",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "SAS",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "EAS",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "UNK",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX.RData"))
  results <- rbind(results,data.frame(dataset = "original",trait = trait,ancestry = "MIX",SL.result))
  
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EUR",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "NonEur",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "AFR",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "SAS",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "EAS",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "UNK",SL.result))
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX_modified.RData"))
  results <- rbind(results,data.frame(dataset = "modified",trait = trait,ancestry = "MIX",SL.result))
}


rm(list=setdiff(ls(), "results"))

library(stringr)

results$method[str_detect(results$method,"CT")] <- "CT"
results$method[str_detect(results$method,"SL_Combined")] <- "SL_Combined"
results$method[str_detect(results$method,"LDPred2")] <- "LDPred2"
results$method[str_detect(results$method,"LASSOSUM2")] <- "LASSOSUM2"
results$method[str_detect(results$method,"CV_plus_RV_STAARO")] <- "CV_plus_RV_STAARO"
results$method[str_detect(results$method,"CV_plus_RV_Burden")] <- "CV_plus_RV_Burden"
results$method[str_detect(results$method,"SlidingWindow_STAARO")] <- "SlidingWindow_STAARO"
results$method[str_detect(results$method,"SlidingWindow_Burden")] <- "SlidingWindow_Burden"
results$method[str_detect(results$method,"GeneCentric_Coding_STAARO")] <- "GeneCentric_Coding_STAARO"
results$method[str_detect(results$method,"GeneCentric_Coding_Burden")] <- "GeneCentric_Coding_Burden"
results$method[str_detect(results$method,"GeneCentric_Noncoding_STAARO")] <- "GeneCentric_Noncoding_STAARO"
results$method[str_detect(results$method,"GeneCentric_Noncoding_Burden")] <- "GeneCentric_Noncoding_Burden"

results$method <- paste0(results$method,"_",results$dataset)

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


ggplot(results[results$ancestry != "EAS",]) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
  #geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("Mod vs OG") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()
