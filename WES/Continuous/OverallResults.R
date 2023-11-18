rm(list = ls())

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

trait <- "BMI"

results <- NULL

for(trait in c("BMI","TC","LDL","HDL","logTG","Height")){
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EUR.RData"))
  ct_EUR <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_NonEur.RData"))
  ct_NonEur <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_AFR.RData"))
  ct_AFR <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_EAS.RData"))
  ct_EAS <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_SAS.RData"))
  ct_SAS <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_MIX.RData"))
  ct_MIX <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/CT/",trait,"_CT_result_UNK.RData"))
  ct_UNK <- ct.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EUR.RData"))
  ldpred2_EUR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_NonEur.RData"))
  ldpred2_NonEur <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_AFR.RData"))
  ldpred2_AFR <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_EAS.RData"))
  ldpred2_EAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_SAS.RData"))
  ldpred2_SAS <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_MIX.RData"))
  ldpred2_MIX <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LDPred2/",trait,"_ldpred2_result_UNK.RData"))
  ldpred2_UNK <- ldpred2.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EUR.RData"))
  LASSOSUM2_EUR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_NonEur.RData"))
  LASSOSUM2_NonEur <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_AFR.RData"))
  LASSOSUM2_AFR <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_EAS.RData"))
  LASSOSUM2_EAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_SAS.RData"))
  LASSOSUM2_SAS <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_MIX.RData"))
  LASSOSUM2_MIX <- LASSOSUM2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/LASSOSUM2/",trait,"_LASSOSUM2_result_UNK.RData"))
  LASSOSUM2_UNK <- LASSOSUM2.result[[1]]
  
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_Eur.RData"))
  sl_EUR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_NonEur.RData"))
  sl_NonEur <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_AFR.RData"))
  sl_AFR <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_EAS.RData"))
  sl_EAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_SAS.RData"))
  sl_SAS <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_MIX.RData"))
  sl_MIX <- SL.result
  
  load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Continuous/Results/Combined_Common_PRS/",trait,"_sl_result_All_UNK.RData"))
  sl_UNK <- SL.result
  
  results_EUR <- rbind(sl_EUR,ct_EUR,LASSOSUM2_EUR,ldpred2_EUR)
  results_EUR$Ancestry <- "EUR"
  results_NonEur <- rbind(sl_NonEur,ct_NonEur,LASSOSUM2_NonEur,ldpred2_NonEur)
  results_NonEur$Ancestry <- "NonEUR"
  results_AFR <- rbind(sl_AFR,ct_AFR,LASSOSUM2_AFR,ldpred2_AFR)
  results_AFR$Ancestry <- "AFR"
  results_EAS <- rbind(sl_EAS,ct_EAS,LASSOSUM2_EAS,ldpred2_EAS)
  results_EAS$Ancestry <- "EAS"
  results_SAS <- rbind(sl_SAS,ct_SAS,LASSOSUM2_SAS,ldpred2_SAS)
  results_SAS$Ancestry <- "SAS"
  results_MIX <- rbind(sl_MIX,ct_MIX,LASSOSUM2_MIX,ldpred2_MIX)
  results_MIX$Ancestry <- "MIX"
  results_UNK <- rbind(sl_UNK,ct_UNK,LASSOSUM2_UNK,ldpred2_UNK)
  results_UNK$Ancestry <- "UNK"
  
  results_tmp <- rbind(results_EUR,results_NonEur,results_AFR,results_EAS,results_SAS,results_MIX,results_UNK)
  results_tmp$Trait <- trait
   
  results <- rbind(results,results_tmp)
}

library(stringr)

results$method[str_detect(results$method,"CT")] <- "CT"
results$method[str_detect(results$method,"SL_Combined")] <- "SL_Combined"
results$method[str_detect(results$method,"LDPred2")] <- "LDPred2"
results$method[str_detect(results$method,"LASSOSUM2")] <- "LASSOSUM2"

library(ggplot2)

ggplot(results) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Common Variants") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()

ggplot(results[results$Ancestry == "EUR",]) +
  geom_bar(aes(x=method, y=r2,fill=method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait)) + 
  ggtitle("Common Variants, EUR") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication() + 
  scale_fill_Publication()
