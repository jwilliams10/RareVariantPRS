rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")

results_70 <- NULL

for(i in 1:length(Y_train)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_Common_PRS/sl_result_All",i,".RData"))
  sl_all <- SL.result
  sl_all$method <- "SL_All"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/CT_result",i,".RData"))
  ct <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/lassosum2_result",i,".RData"))
  lassosum2 <- lassosum2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/ldpred2_result",i,".RData"))
  ldpred2 <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result",i,".RData"))
  GeneCentric_Coding_STAARO <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result",i,".RData"))
  GeneCentric_Coding_Burden <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result",i,".RData"))
  GeneCentric_Noncoding_STAARO <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result",i,".RData"))
  GeneCentric_Noncoding_Burden <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_STAARO_result",i,".RData"))
  SlidingWindow_STAARO <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/SlidingWindow/SlidingWindow_Burden_result",i,".RData"))
  SlidingWindow_Burden <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_STAARO",i,".RData"))
  BestAll_RareVariant_STAARO <- SL.result
  BestAll_RareVariant_STAARO$method <- "SL_All_STAARO"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Combined_RareVariants_PRS/sl_result_All_Burden",i,".RData"))
  BestAll_RareVariant_Burden <- SL.result
  BestAll_RareVariant_Burden$method <- "SL_All_Burden"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/STAARO_All_Result",i,".RData"))
  BestAll_STAARO <- SL.result
  BestAll_STAARO$method <- "CV_RV_STAARO_All"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/Burden_All_Result",i,".RData"))
  BestAll_Burden <- SL.result
  BestAll_Burden$method <- "CV_RV_Burden_All"
  
  
  results_tmp <- rbind(sl_all,ct,lassosum2,ldpred2,
                       GeneCentric_Coding_STAARO,GeneCentric_Coding_Burden,
                       GeneCentric_Noncoding_STAARO,GeneCentric_Noncoding_Burden,
                       SlidingWindow_STAARO,SlidingWindow_Burden,
                       BestAll_RareVariant_STAARO,BestAll_RareVariant_Burden,BestAll_Burden,BestAll_STAARO)
  
  results_70 <- rbind(results_70,results_tmp)
  
  rm(list=setdiff(ls(), c("results_tmp","results_70","i","Y_train"))) 
}

scale <- rep(rep(c("Unscaled","Scaled"),each = 20),3)
causal_prop <- rep(c("Causal Prop. = 0.001","Causal Prop. = 0.01","Causal Prop. = 0.05"),each = 40)

results_70 <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_70$method,R2 = results_70$r2,R2_Low = results_70$r2_low,R2_High = results_70$r2_high)
results_70$Train_Size <- nrow(Y_train[[1]])



rm(list=setdiff(ls(), c("results_70"))) 

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")

results_35 <- NULL

for(i in 1:length(Y_train)){
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_Common_PRS/sl_result_All",i,".RData"))
  sl_all <- SL.result
  sl_all$method <- "SL_All"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/CT_result",i,".RData"))
  ct <- ct.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/lassosum2_result",i,".RData"))
  lassosum2 <- lassosum2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/ldpred2_result",i,".RData"))
  ldpred2 <- ldpred2.result[[1]]
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_STAARO_result",i,".RData"))
  GeneCentric_Coding_STAARO <- r2.result_GeneCentric_Coding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricCoding/GeneCentric_Coding_Burden_result",i,".RData"))
  GeneCentric_Coding_Burden <- r2.result_GeneCentric_Coding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_STAARO_result",i,".RData"))
  GeneCentric_Noncoding_STAARO <- r2.result_GeneCentric_Noncoding_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/GeneCentricNoncoding/GeneCentric_Noncoding_Burden_result",i,".RData"))
  GeneCentric_Noncoding_Burden <- r2.result_GeneCentric_Noncoding_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_STAARO_result",i,".RData"))
  SlidingWindow_STAARO <- r2.result_SlidingWindow_STAARO
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/SlidingWindow/SlidingWindow_Burden_result",i,".RData"))
  SlidingWindow_Burden <- r2.result_SlidingWindow_Burden
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_STAARO",i,".RData"))
  BestAll_RareVariant_STAARO <- SL.result
  BestAll_RareVariant_STAARO$method <- "SL_All_STAARO"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Combined_RareVariants_PRS/sl_result_All_Burden",i,".RData"))
  BestAll_RareVariant_Burden <- SL.result
  BestAll_RareVariant_Burden$method <- "SL_All_Burden"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/STAARO_All_Result",i,".RData"))
  BestAll_STAARO <- SL.result
  BestAll_STAARO$method <- "CV_RV_STAARO_All"
  
  load(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/Burden_All_Result",i,".RData"))
  BestAll_Burden <- SL.result
  BestAll_Burden$method <- "CV_RV_Burden_All"
  
  
  results_tmp <- rbind(sl_all,ct,lassosum2,ldpred2,
                       GeneCentric_Coding_STAARO,GeneCentric_Coding_Burden,
                       GeneCentric_Noncoding_STAARO,GeneCentric_Noncoding_Burden,
                       SlidingWindow_STAARO,SlidingWindow_Burden,
                       BestAll_RareVariant_STAARO,BestAll_RareVariant_Burden,BestAll_Burden,BestAll_STAARO)
  
  results_35 <- rbind(results_35,results_tmp)
  
  rm(list=setdiff(ls(), c("results_tmp","results_35","results_70","i","Y_train"))) 
}

scale <- rep(rep(c("Unscaled","Scaled"),each = 20),3)
causal_prop <- rep(c("Causal Prop. = 0.001","Causal Prop. = 0.01","Causal Prop. = 0.05"),each = 40)

results_35 <- data.frame(Scale = scale, Causal_Prop = causal_prop, Method = results_35$method,R2 = results_35$r2,R2_Low = results_35$r2_low,R2_High = results_35$r2_high)
results_35$Train_Size <- nrow(Y_train[[1]])

results <- rbind(results_35,results_70)
results$Train_Size <- paste0("n = ",results$Train_Size)

rm(list=setdiff(ls(), c("results"))) 

results <- aggregate(.~Method + Scale + Causal_Prop + Train_Size,data = results,mean)

rarevariant_results <- results[results$Method %in% c("GeneCentric_Coding_Burden","GeneCentric_Coding_STAARO","GeneCentric_Noncoding_Burden","GeneCentric_Noncoding_STAARO",
                                                     "SlidingWindow_Burden","SlidingWindow_STAARO","SL_All_STAARO","SL_All_Burden"),]

cv_results <- results[results$Method %in% c("CT","LASSOSUM2","LDPred2","SL_All"),]

overall_results <- results[results$Method %in% c("CT","LASSOSUM2","LDPred2","CV_RV_STAARO_All","CV_RV_Burden_All"),]

####################################################### Plots

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_blank(),
            axis.text = element_text(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

library(ggplot2)

ggplot(rarevariant_results[rarevariant_results$Scale == "Scaled",]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Rare Variant Methods; Scaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication()

ggplot(rarevariant_results[rarevariant_results$Scale == "Unscaled",]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) +
  ggtitle("Rare Variant Methods; Unscaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication()



ggplot(cv_results[cv_results$Scale == "Scaled",]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Common Variant Methods; Scaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication()

ggplot(cv_results[cv_results$Scale == "Unscaled",]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Common Variant Methods; Unscaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication()



ggplot(overall_results[overall_results$Scale == "Scaled",]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Overall Results; Scaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication()

ggplot(overall_results[overall_results$Scale == "Unscaled",]) +
  geom_bar( aes(x=Method, y=R2,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=R2_Low, ymax=R2_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Causal_Prop), vars(Train_Size)) + 
  ggtitle("Overall Results; Unscaled G") + 
  ylab(bquote("R"^"2")) + 
  theme_Publication()
