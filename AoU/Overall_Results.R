rm(list = ls())

library(dplyr)
library(ggpubr)

CTSLEB_Results <- read.csv("/data/williamsjacr/AoU_Results/CTSLEB_Results.csv")
CTSLEB_Results$Method <- "CT-SLEB"
CTSLEB_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CTSLEB_Boot.csv")
CTSLEB_Boot_Results$Method <- "CT-SLEB"
PROSPER_Results <- read.csv("/data/williamsjacr/AoU_Results/PROSPER_Results.csv")
PROSPER_Results$Method <- "PROSPER"
PROSPER_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/PROSPER_Boot.csv")
PROSPER_Boot_Results$Method <- "PROSPER"
JointPRS_Results <- read.csv("/data/williamsjacr/AoU_Results/JointPRS_Results.csv")
JointPRS_Results$Method <- "JointPRS"
JointPRS_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/JointPRS_Boot.csv")
JointPRS_Boot_Results$Method <- "JointPRS"
RICE_CV_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Results.csv")
RICE_CV_Results$Method <- "RICE-CV"
RICE_CV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Boot.csv")
RICE_CV_Boot_Results$Method <- "RICE-CV"
colnames(RICE_CV_Boot_Results) <- colnames(CTSLEB_Boot_Results)
RICE_RV_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Results.csv")
RICE_RV_Results$Method <- "RICE-RV"
RICE_RV_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/RV_Boot.csv")
RICE_RV_Boot_Results$Method <- "RICE-RV"
colnames(RICE_RV_Boot_Results) <- colnames(CTSLEB_Boot_Results)

Comparison_Boot_Results <- read.csv("/data/williamsjacr/AoU_Results/CV_Comparison_Boot.csv")

full_results <- rbind(CTSLEB_Results,PROSPER_Results,JointPRS_Results,RICE_CV_Results,RICE_RV_Results)
full_results_Boot <- rbind(CTSLEB_Boot_Results,PROSPER_Boot_Results,JointPRS_Boot_Results,RICE_CV_Boot_Results,RICE_RV_Boot_Results)

library(ggplot2)
library(boot)

full_results$Method1 <- full_results$Method
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-RV","RICE-CV"))
full_results$Method1[full_results$Method1 == "RICE-RV"] <- "RICE-CV"
full_results$Method1 <- factor(full_results$Method1,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-CV"))

full_results$trait[full_results$trait == "logTG"] <- "log(TG)"
full_results_Boot$trait[full_results_Boot$trait == "logTG"] <- "log(TG)"
Comparison_Boot_Results$trait[Comparison_Boot_Results$trait == "logTG"] <- "log(TG)"
full_results$trait <- factor(full_results$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results_Boot$trait <- factor(full_results_Boot$trait,levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))
full_results$ancestry <- factor(full_results$ancestry,levels = c("AFR","AMR","EAS","EUR","MID","SAS"))

lower_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1,2)] <- paste0(colnames(lower_95)[-c(1,2)],"_Lower")
upper_95 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1,2)] <- paste0(colnames(upper_95)[-c(1,2)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_95)),
                    Method = c(CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method,CI_95$Method),
                    beta_raw_Lower_95 = c(CI_95$beta_raw_AFR_boot_Lower,CI_95$beta_raw_AMR_boot_Lower,CI_95$beta_raw_EAS_boot_Lower,CI_95$beta_raw_EUR_boot_Lower,CI_95$beta_raw_MID_boot_Lower,CI_95$beta_raw_SAS_boot_Lower),
                    beta_raw_Upper_95 = c(CI_95$beta_raw_AFR_boot_Upper,CI_95$beta_raw_AMR_boot_Upper,CI_95$beta_raw_EAS_boot_Upper,CI_95$beta_raw_EUR_boot_Upper,CI_95$beta_raw_MID_boot_Upper,CI_95$beta_raw_SAS_boot_Upper),
                    R2_raw_Lower_95 = c(CI_95$R2_raw_AFR_boot_Lower,CI_95$R2_raw_AMR_boot_Lower,CI_95$R2_raw_EAS_boot_Lower,CI_95$R2_raw_EUR_boot_Lower,CI_95$R2_raw_MID_boot_Lower,CI_95$R2_raw_SAS_boot_Lower),
                    R2_raw_Upper_95 = c(CI_95$R2_raw_AFR_boot_Upper,CI_95$R2_raw_AMR_boot_Upper,CI_95$R2_raw_EAS_boot_Upper,CI_95$R2_raw_EUR_boot_Upper,CI_95$R2_raw_MID_boot_Upper,CI_95$R2_raw_SAS_boot_Upper),
                    beta_adjusted_Lower_95 = c(CI_95$beta_adjusted_AFR_boot_Lower,CI_95$beta_adjusted_AMR_boot_Lower,CI_95$beta_adjusted_EAS_boot_Lower,CI_95$beta_adjusted_EUR_boot_Lower,CI_95$beta_adjusted_MID_boot_Lower,CI_95$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_95 = c(CI_95$beta_adjusted_AFR_boot_Upper,CI_95$beta_adjusted_AMR_boot_Upper,CI_95$beta_adjusted_EAS_boot_Upper,CI_95$beta_adjusted_EUR_boot_Upper,CI_95$beta_adjusted_MID_boot_Upper,CI_95$beta_adjusted_SAS_boot_Upper),
                    R2_adjusted_Lower_95 = c(CI_95$R2_adjusted_AFR_boot_Lower,CI_95$R2_adjusted_AMR_boot_Lower,CI_95$R2_adjusted_EAS_boot_Lower,CI_95$R2_adjusted_EUR_boot_Lower,CI_95$R2_adjusted_MID_boot_Lower,CI_95$R2_adjusted_SAS_boot_Lower),
                    R2_adjusted_Upper_95 = c(CI_95$R2_adjusted_AFR_boot_Upper,CI_95$R2_adjusted_AMR_boot_Upper,CI_95$R2_adjusted_EAS_boot_Upper,CI_95$R2_adjusted_EUR_boot_Upper,CI_95$R2_adjusted_MID_boot_Upper,CI_95$R2_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_95)

lower_95 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.025)})
colnames(lower_95)[-c(1)] <- paste0(colnames(lower_95)[-c(1)],"_Lower")
upper_95 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.975)})
colnames(upper_95)[-c(1)] <- paste0(colnames(upper_95)[-c(1)],"_Upper")
CI_95 <- inner_join(lower_95,upper_95)
CI_95 <- data.frame(trait = c(CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait,CI_95$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_95)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CTSLEB_Lower_95 = c(CI_95$R2_raw_AFR_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_AMR_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_EAS_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_EUR_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_MID_RICE_vs_CTSLEB_Lower,CI_95$R2_raw_SAS_RICE_vs_CTSLEB_Lower),
                    R2_raw_RICE_vs_CTSLEB_Upper_95 = c(CI_95$R2_raw_AFR_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_AMR_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_EAS_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_EUR_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_MID_RICE_vs_CTSLEB_Upper,CI_95$R2_raw_SAS_RICE_vs_CTSLEB_Upper),
                    R2_raw_RICE_vs_PROSPER_Lower_95 = c(CI_95$R2_raw_AFR_RICE_vs_PROSPER_Lower,CI_95$R2_raw_AMR_RICE_vs_PROSPER_Lower,CI_95$R2_raw_EAS_RICE_vs_PROSPER_Lower,CI_95$R2_raw_EUR_RICE_vs_PROSPER_Lower,CI_95$R2_raw_MID_RICE_vs_PROSPER_Lower,CI_95$R2_raw_SAS_RICE_vs_PROSPER_Lower),
                    R2_raw_RICE_vs_PROSPER_Upper_95 = c(CI_95$R2_raw_AFR_RICE_vs_PROSPER_Upper,CI_95$R2_raw_AMR_RICE_vs_PROSPER_Upper,CI_95$R2_raw_EAS_RICE_vs_PROSPER_Upper,CI_95$R2_raw_EUR_RICE_vs_PROSPER_Upper,CI_95$R2_raw_MID_RICE_vs_PROSPER_Upper,CI_95$R2_raw_SAS_RICE_vs_PROSPER_Upper),
                    R2_raw_RICE_vs_JointPRS_Lower_95 = c(CI_95$R2_raw_AFR_RICE_vs_JointPRS_Lower,CI_95$R2_raw_AMR_RICE_vs_JointPRS_Lower,CI_95$R2_raw_EAS_RICE_vs_JointPRS_Lower,CI_95$R2_raw_EUR_RICE_vs_JointPRS_Lower,CI_95$R2_raw_MID_RICE_vs_JointPRS_Lower,CI_95$R2_raw_SAS_RICE_vs_JointPRS_Lower),
                    R2_raw_RICE_vs_JointPRS_Upper_95 = c(CI_95$R2_raw_AFR_RICE_vs_JointPRS_Upper,CI_95$R2_raw_AMR_RICE_vs_JointPRS_Upper,CI_95$R2_raw_EAS_RICE_vs_JointPRS_Upper,CI_95$R2_raw_EUR_RICE_vs_JointPRS_Upper,CI_95$R2_raw_MID_RICE_vs_JointPRS_Upper,CI_95$R2_raw_SAS_RICE_vs_JointPRS_Upper),
                    R2_adjusted_RICE_vs_CTSLEB_Lower_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_AMR_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_EAS_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_EUR_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_MID_RICE_vs_CTSLEB_Lower,CI_95$R2_adjusted_SAS_RICE_vs_CTSLEB_Lower),
                    R2_adjusted_RICE_vs_CTSLEB_Upper_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_AMR_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_EAS_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_EUR_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_MID_RICE_vs_CTSLEB_Upper,CI_95$R2_adjusted_SAS_RICE_vs_CTSLEB_Upper),
                    R2_adjusted_RICE_vs_PROSPER_Lower_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_AMR_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_EAS_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_EUR_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_MID_RICE_vs_PROSPER_Lower,CI_95$R2_adjusted_SAS_RICE_vs_PROSPER_Lower),
                    R2_adjusted_RICE_vs_PROSPER_Upper_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_AMR_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_EAS_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_EUR_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_MID_RICE_vs_PROSPER_Upper,CI_95$R2_adjusted_SAS_RICE_vs_PROSPER_Upper),
                    R2_adjusted_RICE_vs_JointPRS_Lower_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_AMR_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_EAS_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_EUR_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_MID_RICE_vs_JointPRS_Lower,CI_95$R2_adjusted_SAS_RICE_vs_JointPRS_Lower),
                    R2_adjusted_RICE_vs_JointPRS_Upper_95 = c(CI_95$R2_adjusted_AFR_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_AMR_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_EAS_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_EUR_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_MID_RICE_vs_JointPRS_Upper,CI_95$R2_adjusted_SAS_RICE_vs_JointPRS_Upper))
full_results <- left_join(full_results,CI_95)

lower_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.005)})
colnames(lower_99)[-c(1,2)] <- paste0(colnames(lower_99)[-c(1,2)],"_Lower")
upper_99 <- aggregate(.~trait + Method,data = full_results_Boot,function(x){quantile(x,0.995)})
colnames(upper_99)[-c(1,2)] <- paste0(colnames(upper_99)[-c(1,2)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_99)),
                    Method = c(CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method,CI_99$Method),
                    beta_raw_Lower_99 = c(CI_99$beta_raw_AFR_boot_Lower,CI_99$beta_raw_AMR_boot_Lower,CI_99$beta_raw_EAS_boot_Lower,CI_99$beta_raw_EUR_boot_Lower,CI_99$beta_raw_MID_boot_Lower,CI_99$beta_raw_SAS_boot_Lower),
                    beta_raw_Upper_99 = c(CI_99$beta_raw_AFR_boot_Upper,CI_99$beta_raw_AMR_boot_Upper,CI_99$beta_raw_EAS_boot_Upper,CI_99$beta_raw_EUR_boot_Upper,CI_99$beta_raw_MID_boot_Upper,CI_99$beta_raw_SAS_boot_Upper),
                    R2_raw_Lower_99 = c(CI_99$R2_raw_AFR_boot_Lower,CI_99$R2_raw_AMR_boot_Lower,CI_99$R2_raw_EAS_boot_Lower,CI_99$R2_raw_EUR_boot_Lower,CI_99$R2_raw_MID_boot_Lower,CI_99$R2_raw_SAS_boot_Lower),
                    R2_raw_Upper_99 = c(CI_99$R2_raw_AFR_boot_Upper,CI_99$R2_raw_AMR_boot_Upper,CI_99$R2_raw_EAS_boot_Upper,CI_99$R2_raw_EUR_boot_Upper,CI_99$R2_raw_MID_boot_Upper,CI_99$R2_raw_SAS_boot_Upper),
                    beta_adjusted_Lower_99 = c(CI_99$beta_adjusted_AFR_boot_Lower,CI_99$beta_adjusted_AMR_boot_Lower,CI_99$beta_adjusted_EAS_boot_Lower,CI_99$beta_adjusted_EUR_boot_Lower,CI_99$beta_adjusted_MID_boot_Lower,CI_99$beta_adjusted_SAS_boot_Lower),
                    beta_adjusted_Upper_99 = c(CI_99$beta_adjusted_AFR_boot_Upper,CI_99$beta_adjusted_AMR_boot_Upper,CI_99$beta_adjusted_EAS_boot_Upper,CI_99$beta_adjusted_EUR_boot_Upper,CI_99$beta_adjusted_MID_boot_Upper,CI_99$beta_adjusted_SAS_boot_Upper),
                    R2_adjusted_Lower_99 = c(CI_99$R2_adjusted_AFR_boot_Lower,CI_99$R2_adjusted_AMR_boot_Lower,CI_99$R2_adjusted_EAS_boot_Lower,CI_99$R2_adjusted_EUR_boot_Lower,CI_99$R2_adjusted_MID_boot_Lower,CI_99$R2_adjusted_SAS_boot_Lower),
                    R2_adjusted_Upper_99 = c(CI_99$R2_adjusted_AFR_boot_Upper,CI_99$R2_adjusted_AMR_boot_Upper,CI_99$R2_adjusted_EAS_boot_Upper,CI_99$R2_adjusted_EUR_boot_Upper,CI_99$R2_adjusted_MID_boot_Upper,CI_99$R2_adjusted_SAS_boot_Upper)) 
full_results <- left_join(full_results,CI_99)

lower_99 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.025)})
colnames(lower_99)[-c(1)] <- paste0(colnames(lower_99)[-c(1)],"_Lower")
upper_99 <- aggregate(.~trait,data = Comparison_Boot_Results,function(x){quantile(x,0.975)})
colnames(upper_99)[-c(1)] <- paste0(colnames(upper_99)[-c(1)],"_Upper")
CI_99 <- inner_join(lower_99,upper_99)
CI_99 <- data.frame(trait = c(CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait,CI_99$trait),
                    ancestry = rep(c("AFR","AMR","EAS","EUR","MID","SAS"),each = nrow(CI_99)),
                    Method = "RICE-CV",
                    R2_raw_RICE_vs_CTSLEB_Lower_99 = c(CI_99$R2_raw_AFR_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_AMR_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_EAS_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_EUR_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_MID_RICE_vs_CTSLEB_Lower,CI_99$R2_raw_SAS_RICE_vs_CTSLEB_Lower),
                    R2_raw_RICE_vs_CTSLEB_Upper_99 = c(CI_99$R2_raw_AFR_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_AMR_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_EAS_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_EUR_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_MID_RICE_vs_CTSLEB_Upper,CI_99$R2_raw_SAS_RICE_vs_CTSLEB_Upper),
                    R2_raw_RICE_vs_PROSPER_Lower_99 = c(CI_99$R2_raw_AFR_RICE_vs_PROSPER_Lower,CI_99$R2_raw_AMR_RICE_vs_PROSPER_Lower,CI_99$R2_raw_EAS_RICE_vs_PROSPER_Lower,CI_99$R2_raw_EUR_RICE_vs_PROSPER_Lower,CI_99$R2_raw_MID_RICE_vs_PROSPER_Lower,CI_99$R2_raw_SAS_RICE_vs_PROSPER_Lower),
                    R2_raw_RICE_vs_PROSPER_Upper_99 = c(CI_99$R2_raw_AFR_RICE_vs_PROSPER_Upper,CI_99$R2_raw_AMR_RICE_vs_PROSPER_Upper,CI_99$R2_raw_EAS_RICE_vs_PROSPER_Upper,CI_99$R2_raw_EUR_RICE_vs_PROSPER_Upper,CI_99$R2_raw_MID_RICE_vs_PROSPER_Upper,CI_99$R2_raw_SAS_RICE_vs_PROSPER_Upper),
                    R2_raw_RICE_vs_JointPRS_Lower_99 = c(CI_99$R2_raw_AFR_RICE_vs_JointPRS_Lower,CI_99$R2_raw_AMR_RICE_vs_JointPRS_Lower,CI_99$R2_raw_EAS_RICE_vs_JointPRS_Lower,CI_99$R2_raw_EUR_RICE_vs_JointPRS_Lower,CI_99$R2_raw_MID_RICE_vs_JointPRS_Lower,CI_99$R2_raw_SAS_RICE_vs_JointPRS_Lower),
                    R2_raw_RICE_vs_JointPRS_Upper_99 = c(CI_99$R2_raw_AFR_RICE_vs_JointPRS_Upper,CI_99$R2_raw_AMR_RICE_vs_JointPRS_Upper,CI_99$R2_raw_EAS_RICE_vs_JointPRS_Upper,CI_99$R2_raw_EUR_RICE_vs_JointPRS_Upper,CI_99$R2_raw_MID_RICE_vs_JointPRS_Upper,CI_99$R2_raw_SAS_RICE_vs_JointPRS_Upper),
                    R2_adjusted_RICE_vs_CTSLEB_Lower_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_AMR_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_EAS_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_EUR_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_MID_RICE_vs_CTSLEB_Lower,CI_99$R2_adjusted_SAS_RICE_vs_CTSLEB_Lower),
                    R2_adjusted_RICE_vs_CTSLEB_Upper_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_AMR_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_EAS_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_EUR_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_MID_RICE_vs_CTSLEB_Upper,CI_99$R2_adjusted_SAS_RICE_vs_CTSLEB_Upper),
                    R2_adjusted_RICE_vs_PROSPER_Lower_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_AMR_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_EAS_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_EUR_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_MID_RICE_vs_PROSPER_Lower,CI_99$R2_adjusted_SAS_RICE_vs_PROSPER_Lower),
                    R2_adjusted_RICE_vs_PROSPER_Upper_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_AMR_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_EAS_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_EUR_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_MID_RICE_vs_PROSPER_Upper,CI_99$R2_adjusted_SAS_RICE_vs_PROSPER_Upper),
                    R2_adjusted_RICE_vs_JointPRS_Lower_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_AMR_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_EAS_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_EUR_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_MID_RICE_vs_JointPRS_Lower,CI_99$R2_adjusted_SAS_RICE_vs_JointPRS_Lower),
                    R2_adjusted_RICE_vs_JointPRS_Upper_99 = c(CI_99$R2_adjusted_AFR_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_AMR_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_EAS_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_EUR_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_MID_RICE_vs_JointPRS_Upper,CI_99$R2_adjusted_SAS_RICE_vs_JointPRS_Upper))
full_results <- left_join(full_results,CI_99)

full_results_stacked <- rbind(data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_raw, lower_95 = full_results$beta_raw_Lower_95, upper_95 = full_results$beta_raw_Upper_95,method = full_results$Method,Standardization = "Within Genetically-Inferred Ancestries"),
                              data.frame(trait = full_results$trait, ancestry = full_results$ancestry,beta = full_results$beta_adjusted, lower_95 = full_results$beta_adjusted_Lower_95, upper_95 = full_results$beta_adjusted_Upper_95,method = full_results$Method,Standardization = "Using PCs 1-5"))


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

g1 <- ggplot(data=full_results_stacked[full_results_stacked$trait %in% c("BMI","Height"),], aes(x=method, y=beta, ymin=lower_95, ymax=upper_95,color = Standardization)) +
  geom_pointrange(position=position_dodge(width=.25),size = 0.2) + 
  facet_grid(vars(trait), vars(ancestry)) + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Beta of PRS per SD") +
  theme_Publication() + 
  scale_fill_Publication()
ggsave(paste0("AoU_WES_Continuous_","A","_Raw_vs_AncestryAdjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)

g1 <- ggplot(data=full_results_stacked[full_results_stacked$trait %in% c("HDL","LDL"),], aes(x=method, y=beta, ymin=lower_95, ymax=upper_95,color = Standardization)) +
  geom_pointrange(position=position_dodge(width=.25),size = 0.2) + 
  facet_grid(vars(trait), vars(ancestry)) + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Beta of PRS per SD") +
  theme_Publication() + 
  scale_fill_Publication()
ggsave(paste0("AoU_WES_Continuous_","B","_Raw_vs_AncestryAdjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)


g1 <- ggplot(data=full_results_stacked[full_results_stacked$trait %in% c("log(TG)","TC"),], aes(x=method, y=beta, ymin=lower_95, ymax=upper_95,color = Standardization)) +
  geom_pointrange(position=position_dodge(width=.25),size = 0.2) + 
  facet_grid(vars(trait), vars(ancestry)) + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Method") + ylab("Beta of PRS per SD") +
  theme_Publication() + 
  scale_fill_Publication()
ggsave(paste0("AoU_WES_Continuous_","C","_Raw_vs_AncestryAdjusted_Beta.png"),g1,width=10, height=6.18047,dpi = 300)

# ggplot(full_results) +
#   geom_bar(aes(x=Method, y=abs(beta_raw),fill=Method), stat="identity", alpha=0.7) +
#   facet_grid(vars(trait), vars(ancestry)) + 
#   ggtitle("AoU WGS Raw PRS Results") + 
#   ylab("Beta of PRS per SD") + 
#   ylim(0,0.6) +
#   theme_Publication() + 
#   scale_fill_Publication()

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

full_results$beta_adjusted[full_results$beta_adjusted < 0] <- 0
full_results$beta_raw[full_results$beta_raw < 0] <- 0

full_results$group1 <- "RICE-CV"
full_results$group2 <- "RICE-CV"
full_results$p.signif_beta <- ""
full_results$p.signif_beta[full_results$Method == "RICE-CV"] <- ifelse(full_results$beta_adjusted_Lower_99[full_results$Method == "RICE-RV"] > 0,"***",ifelse(full_results$beta_adjusted_Lower_95[full_results$Method == "RICE-RV"] > 0,"**",""))
full_results$position <- NA
full_results$position[full_results$Method == "RICE-CV"] <- full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"] + 0.03
ylim <- max(c(full_results$beta_adjusted[full_results$Method == "RICE-CV"] + full_results$beta_adjusted[full_results$Method == "RICE-RV"],full_results$beta_adjusted)) + 0.05

full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE-RV","RICE-CV"))

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method1, y=abs(beta_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +
  facet_grid(vars(trait), vars(ancestry)) +
  ggtitle("AoU PRS results for Continuous Traits Across Six Ancestral Groups") +
  ylab("Beta of PRS per SD") +
  ylim(0,ylim) +
  stat_pvalue_manual(full_results,
                     label = "p.signif_beta",
                     y.position = "position",
                     size = 2) +
  theme_Publication() +
  scale_fill_Publication()

ggsave(paste0("AoU_WES_Continuous_Adjusted_Beta.png"),g2,width=10, height=6.18047,dpi = 300)

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#5EBD3E","#FFB900","#F78200","#973999","#009cdf")), ...)
}

full_results$Method <- as.character(full_results$Method)
full_results <- full_results[full_results$Method != "RICE-RV",]
full_results$Method[full_results$Method == "RICE-CV"] <- "RICE"
full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE"))

full_results$group1 <- "RICE"
full_results$group2 <- "RICE"
full_results$p.signif_beta1 <- ""
full_results$p.signif_beta2 <- ""

for(trait in c("BMI","Height","LDL","log(TG)","TC","HDL")){
  for(anc in c("AFR","EUR","SAS","AMR","MID","EAS")){
    tmp <- full_results[full_results$ancestry == anc & full_results$trait == trait,]
    max_R2_notRICE <- max(tmp$R2_adjusted[tmp$Method != "RICE"])
    Best_Method <- tmp$Method[tmp$R2_adjusted == max_R2_notRICE]
    Improved_R2 <- round((tmp$R2_adjusted[tmp$Method == "RICE"]/max_R2_notRICE - 1)*100,digits = 2)
    
    if(Best_Method == "CT-SLEB"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_CTSLEB_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else if(Best_Method == "JointPRS"){
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_JointPRS_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }else{
      full_results$p.signif_beta1[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_99[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_95[tmp$Method == "RICE"] > 0,paste0(Improved_R2,"%"),""))
      full_results$p.signif_beta2[full_results$ancestry == anc & full_results$trait == trait & full_results$Method == "RICE"] <- ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_99[tmp$Method == "RICE"] > 0,paste0("***"),ifelse(tmp$R2_adjusted_RICE_vs_PROSPER_Lower_95[tmp$Method == "RICE"] > 0,paste0("**"),""))
    }
  }
}

full_results$position1 <- NA
full_results$position2 <- NA
full_results$position1[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.01
full_results$position2[full_results$Method == "RICE"] <- full_results$R2_adjusted[full_results$Method == "RICE"] + 0.07
ylim <- max(c(full_results$R2_adjusted)) + 0.08

full_results$Method <- factor(full_results$Method,levels = c("CT-SLEB","JointPRS","PROSPER","RICE"))

g2 <- ggplot(full_results) +
  geom_bar(aes(x=Method, y=abs(R2_adjusted),fill=Method), stat="identity", alpha=0.7) +
  # geom_errorbar( aes(x=Method, ymin=r2_low, ymax=r2_high), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(trait), vars(ancestry)) + 
  ggtitle("AoU PRS results for Continuous Traits Across Six Ancestral Groups") + 
  ylab("R2") + 
  ylim(0,ylim) +
  stat_pvalue_manual(full_results,
                     label = "p.signif_beta1",
                     y.position = "position1",
                     size = 2) +
  stat_pvalue_manual(full_results,
                     label = "p.signif_beta2",
                     y.position = "position2",
                     size = 2) +
  theme_Publication() + 
  scale_fill_Publication()

# ggsave(paste0("Desktop/RareVariantPRS_Results/Figures/UKB_Imputed_Continuous_Adjusted_R2.png"),g2,width=10, height=6.18047,dpi = 300)
ggsave(paste0("AoU_WES_Continuous_Adjusted_R2.png"),g2,width=10, height=6.18047,dpi = 300)
