rm(list = ls())

library(dplyr)

index_mat <- NULL

causalprop_vec <- c(0.2,0.05,0.01,0.001,0.0005)
scale <- c(0,1)

count <- 1

for(j in 1:length(causalprop_vec)){
  for(q in 1:length(scale)){
    for(l in 1:100){
      index_mat <- rbind(index_mat,data.frame(i = count,Causal_Prop = causalprop_vec[j],Scale = scale[q]))
      count <- count + 1
    }
  }
}

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_70 <- NULL

for(i in 1:length(Y_train)){

  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"

  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"

  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"

  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation1/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"

  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)

  results_70 <- rbind(results_70,betas_tmp)

  rm(list=setdiff(ls(), c("results_70","i","Y_train","index_mat")))
}

results_70 <- inner_join(results_70,index_mat)
results_70$Causal_Prop <- as.character(results_70$Causal_Prop)
results_70$Causal_Prop[results_70$Causal_Prop == "5e-04"] <- "0.0005"
results_70$Causal_Prop <- paste0("Causal Prop. ",results_70$Causal_Prop)

results_70$Scale <- as.character(results_70$Scale)
results_70$Scale[results_70$Scale == "0"] <- "Unscaled"
results_70$Scale[results_70$Scale == "1"] <- "Scaled"

results_70 <- data.frame(Scale = results_70$Scale, Causal_Prop = results_70$Causal_Prop, Method = results_70$Method,Ancestry = results_70$ancestry,
                         Beta = results_70$beta_adjusted,SE_Beta = results_70$beta_se_adjusted,Lower_Beta = results_70$beta_lower_adjusted,Upper_Beta_Beta = results_70$beta_upper_adjusted,
                         R2 = results_70$R2_adjusted,SE_Beta = results_70$R2_se_adjusted,Lower_Beta = results_70$R2_lower_adjusted,Upper_Beta_Beta = results_70$R2_upper_adjusted)
results_70$Train_Size <- nrow(Y_train[[1]])







load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_35 <- NULL

for(i in 1:length(Y_train)){
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation2/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  results_35 <- rbind(results_35,betas_tmp)

  rm(list=setdiff(ls(), c("results_70","results_35","i","Y_train","index_mat")))
}

results_35 <- inner_join(results_35,index_mat)
results_35$Causal_Prop <- as.character(results_35$Causal_Prop)
results_35$Causal_Prop[results_35$Causal_Prop == "5e-04"] <- "0.0005"
results_35$Causal_Prop <- paste0("Causal Prop. ",results_35$Causal_Prop)

results_35$Scale <- as.character(results_35$Scale)
results_35$Scale[results_35$Scale == "0"] <- "Unscaled"
results_35$Scale[results_35$Scale == "1"] <- "Scaled"

results_35 <- data.frame(Scale = results_35$Scale, Causal_Prop = results_35$Causal_Prop, Method = results_35$Method,Ancestry = results_35$ancestry,
                        Beta = results_35$beta_adjusted,SE_Beta = results_35$beta_se_adjusted,Lower_Beta = results_35$beta_lower_adjusted,Upper_Beta_Beta = results_35$beta_upper_adjusted,
                        R2 = results_35$R2_adjusted,SE_Beta = results_35$R2_se_adjusted,Lower_Beta = results_35$R2_lower_adjusted,Upper_Beta_Beta = results_35$R2_upper_adjusted)
results_35$Train_Size <- nrow(Y_train[[1]])





load("/data/williamsjacr/UKB_WES_Simulation/Simulation3/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_rareprop_70 <- NULL

for(i in 1:length(Y_train)){
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation3/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  
  results_rareprop_70 <- rbind(results_rareprop_70,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_rareprop_70","results_70","results_35","i","Y_train","index_mat")))
}

results_rareprop_70 <- inner_join(results_rareprop_70,index_mat)
results_rareprop_70$Causal_Prop <- as.character(results_rareprop_70$Causal_Prop)
results_rareprop_70$Causal_Prop[results_rareprop_70$Causal_Prop == "5e-04"] <- "0.0005"
results_rareprop_70$Causal_Prop <- paste0("Causal Prop. ",results_rareprop_70$Causal_Prop)

results_rareprop_70$Scale <- as.character(results_rareprop_70$Scale)
results_rareprop_70$Scale[results_rareprop_70$Scale == "0"] <- "Unscaled"
results_rareprop_70$Scale[results_rareprop_70$Scale == "1"] <- "Scaled"

results_rareprop_70 <- data.frame(Scale = results_rareprop_70$Scale, Causal_Prop = results_rareprop_70$Causal_Prop, Method = results_rareprop_70$Method,Ancestry = results_rareprop_70$ancestry,
                                  Beta = results_rareprop_70$beta_adjusted,SE_Beta = results_rareprop_70$beta_se_adjusted,Lower_Beta = results_rareprop_70$beta_lower_adjusted,Upper_Beta_Beta = results_rareprop_70$beta_upper_adjusted,
                                  R2 = results_rareprop_70$R2_adjusted,SE_Beta = results_rareprop_70$R2_se_adjusted,Lower_Beta = results_rareprop_70$R2_lower_adjusted,Upper_Beta_Beta = results_rareprop_70$R2_upper_adjusted)
results_rareprop_70$Train_Size <- nrow(Y_train[[1]])







load("/data/williamsjacr/UKB_WES_Simulation/Simulation4/simulated_data/phenotypes/Y_Train.RData")

i <- 1

results_rareprop_35 <- NULL

for(i in 1:length(Y_train)){
  Best_Betas_CT <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/CT/Best_Betas",i,".csv"))
  Best_Betas_CT$Method <- "CT"
  
  Best_Betas_LDPred <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/LDPred2/Best_Betas",i,".csv"))
  Best_Betas_LDPred$Method <- "LDPred"
  
  Best_Betas_LASSOSum <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/LASSOSUM2/Best_Betas",i,".csv"))
  Best_Betas_LASSOSum$Method <- "LASSOSum"
  
  Best_Betas_RICECV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/Common_plus_RareVariants/CV_Best_Betas",i,".csv"))
  Best_Betas_RICECV$Method <- "RICE-CV"
  
  Best_Betas_RICERV <- read.csv(paste0("/data/williamsjacr/UKB_WES_Simulation/Simulation4/Results/Common_plus_RareVariants/RV_Best_Betas",i,".csv"))
  Best_Betas_RICERV$Method <- "RICE-RV"
  
  betas_tmp <- rbind(Best_Betas_CT,Best_Betas_LDPred,Best_Betas_LASSOSum,Best_Betas_RICECV,Best_Betas_RICERV)
  results_rareprop_35 <- rbind(results_rareprop_35,betas_tmp)
  
  rm(list=setdiff(ls(), c("results_rareprop_35","results_rareprop_70","results_70","results_35","i","Y_train","index_mat")))
}

results_rareprop_35 <- inner_join(results_rareprop_35,index_mat)
results_rareprop_35$Causal_Prop <- as.character(results_rareprop_35$Causal_Prop)
results_rareprop_35$Causal_Prop[results_rareprop_35$Causal_Prop == "5e-04"] <- "0.0005"
results_rareprop_35$Causal_Prop <- paste0("Causal Prop. ",results_rareprop_35$Causal_Prop)

results_rareprop_35$Scale <- as.character(results_rareprop_35$Scale)
results_rareprop_35$Scale[results_rareprop_35$Scale == "0"] <- "Unscaled"
results_rareprop_35$Scale[results_rareprop_35$Scale == "1"] <- "Scaled"

results_rareprop_35 <- data.frame(Scale = results_rareprop_35$Scale, Causal_Prop = results_rareprop_35$Causal_Prop, Method = results_rareprop_35$Method,Ancestry = results_rareprop_35$ancestry,
                                  Beta = results_rareprop_35$beta_adjusted,SE_Beta = results_rareprop_35$beta_se_adjusted,Lower_Beta = results_rareprop_35$beta_lower_adjusted,Upper_Beta_Beta = results_rareprop_35$beta_upper_adjusted,
                                  R2 = results_rareprop_35$R2_adjusted,SE_Beta = results_rareprop_35$R2_se_adjusted,Lower_Beta = results_rareprop_35$R2_lower_adjusted,Upper_Beta_Beta = results_rareprop_35$R2_upper_adjusted)
results_rareprop_35$Train_Size <- nrow(Y_train[[1]])







results <- rbind(results_35,results_70)
results_rareprop <- rbind(results_rareprop_35,results_rareprop_70)

results$Train_Size <- format(results$Train_Size,big.mark=",", trim=TRUE)
results_rareprop$Train_Size <- format(results_rareprop$Train_Size,big.mark=",", trim=TRUE)

results$Train_Size <- paste0("n = ",results$Train_Size)
results_rareprop$Train_Size <- paste0("n = ",results_rareprop$Train_Size)

rm(list=setdiff(ls(), c("results","results_rareprop")))

results$Beta[results$Beta < 0 & results$Method %in% c("LDPred","LASSOSum")] <- -1*results$Beta[results$Beta < 0 & results$Method %in% c("LDPred","LASSOSum")]
results$Beta[results$Beta < 0] <- 0

results_rareprop$Beta[results_rareprop$Beta < 0 & results_rareprop$Method %in% c("LDPred","LASSOSum")] <- -1*results_rareprop$Beta[results_rareprop$Beta < 0 & results_rareprop$Method %in% c("LDPred","LASSOSum")]
results_rareprop$Beta[results_rareprop$Beta < 0] <- 0

results <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results,mean)
results_rareprop <- aggregate(.~Method + Scale + Causal_Prop + Train_Size + Ancestry,data = results_rareprop,mean)

overall_results <- results[results$Method %in% c("CT","LASSOSum","LDPred","RICE-CV","RICE-RV"),]
overall_results_rareprop <- results_rareprop[results_rareprop$Method %in% c("CT","LASSOSum","LDPred","RICE-CV","RICE-RV"),]

# overall_results <- read.csv("Desktop/RareVariantPRS_Results/Overall_Results_SimStudy.csv")

overall_results <- overall_results[overall_results$Ancestry %in% c("AFR","EUR","SAS","AMR"),]
overall_results$Method[overall_results$Method == "RICE-CV"] <- "RICE-CV" 
overall_results$Method[overall_results$Method == "RICE-RV"] <- "RICE-RV" 
overall_results$Method[overall_results$Method == "LDPred"] <- "LDpred2"
overall_results$Method[overall_results$Method == "LASSOSum"] <- "Lassosum2"

overall_results$Method1 <- overall_results$Method
overall_results$Method <- factor(overall_results$Method,levels = c("CT","Lassosum2","LDpred2","RICE-RV","RICE-CV"))
overall_results$Method1[overall_results$Method1 == "RICE-RV"] <- "RICE-CV"
overall_results$Method1 <- factor(overall_results$Method1,levels = c("CT","Lassosum2","LDpred2","RICE-CV"))
overall_results$Ancestry <- factor(overall_results$Ancestry,levels = c("AFR","AMR","EUR","SAS"))


####################################################### Plots

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

library(ggplot2)
library(dplyr)

g1 <- ggplot(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  ylim(0,0.4)+
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Scaled_49173_Adjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)

g2 <- ggplot(overall_results[(overall_results$Scale == "Scaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  ylim(0,0.4) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Scaled_98343_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)

g3 <- ggplot(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 49,173") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 49,173") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  ylim(0,0.4) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Unscaled_49173_Adjusted_Beta.png"),g3,width=10, height=6.18047,dpi = 300)

g4 <- ggplot(overall_results[(overall_results$Scale == "Unscaled") & (overall_results$Train_Size == "n = 98,343") & (overall_results$Causal_Prop %in% c("Causal Prop. 0.2","Causal Prop. 0.05","Causal Prop. 0.01")),]) +
  geom_bar(aes(x=Method1, y=abs(Beta),fill=Method), stat="identity", alpha=0.7) +
  facet_grid(vars(Causal_Prop), vars(Ancestry)) + 
  ggtitle("Simulation using UKB WES with Training Sample Size of 98,343") + 
  ylab("Beta of PRS per SD") + 
  theme_Publication() + 
  ylim(0,0.4) +
  scale_fill_Publication()

ggsave(paste0("UKB_Simulation_Unscaled_98343_Adjusted_Beta.png"),g4,width=10, height=6.18047,dpi = 300)

