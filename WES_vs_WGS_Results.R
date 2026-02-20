# =============================================================================
# [RICE-ANNOTATION] WES_vs_WGS_Results.R
# Purpose: Compares UKB Imputed+WES vs UKB WGS PRS performance (Manuscript Figure 6; Supplementary Figures 15–16).
#
# Paper linkage:
#   - Manuscript: Figure/table generation and cross-platform comparisons.
#   - Supplementary Figures/Data: see script-specific outputs and top-level README.
#
# Notes:
#   - Annotations are intended to point readers to Manuscript, Supplementary Data, and Supplementary Figures.
#   - Many scripts contain environment-specific paths (HPC/DNAnexus/AoU workbench). Update paths as needed for your setup.
#   - See README.md in this directory for expected inputs/outputs and run order.
# =============================================================================
## ============================================================
## UKB PRS Plots (DROP pure WES from plots)
## Keeps: "Imputed + WES" vs "WGS"
## Methods: RICE-CV and RICE-RV
## Outputs: EUR (continuous + binary), AFR/AMR/SAS (continuous)
## ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(cowplot)
library(RColorBrewer)

traits_cont <- c("BMI","TC","HDL","LDL","logTG","Height")
traits_bin  <- c("Asthma","Breast","CAD","Prostate","T2D")

## ------------------------------------------------------------
## Read Imputed + WES results (Continuous)
## ------------------------------------------------------------
Imputed_Results_Continuous <- NULL
for(trait in traits_cont){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  Imputed_Results_Continuous <- rbind(Imputed_Results_Continuous,
                                      rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
Imputed_Results_Continuous$Data_Type <- "Imputed + WES"

## ------------------------------------------------------------
## Read Imputed + WES results (Binary)
## ------------------------------------------------------------
Imputed_Results_Binary <- NULL
for(trait in traits_bin){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LDpred2/",trait,"Best_Betas.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/LASSOsum2/",trait,"Best_Betas.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Imputed/Results/Common_plus_RareVariants/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  Imputed_Results_Binary <- rbind(Imputed_Results_Binary,
                                  rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
Imputed_Results_Binary$Data_Type <- "Imputed + WES"

## ------------------------------------------------------------
## Read WGS results (Binary)
## ------------------------------------------------------------
WGS_Results_Binary <- NULL
for(trait in traits_bin){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results_Binary/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  WGS_Results_Binary <- rbind(WGS_Results_Binary,
                              rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
WGS_Results_Binary$Data_Type <- "WGS"

## ------------------------------------------------------------
## Read WGS results (Continuous)
## ------------------------------------------------------------
WGS_Results_Continuous <- NULL
for(trait in traits_cont){
  CT_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/CT/",trait,"Best_Betas.csv"))
  CT_Results$Method <- "CT"
  
  LDPred2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  LDPred2_Results$Method <- "LDpred2"
  
  LASSOSUM2_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  LASSOSUM2_Results$Method <- "Lassosum2"
  
  RICE_CV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/CV_",trait,"Best_Betas.csv"))
  RICE_CV_Results$Method <- "RICE-CV"
  
  RICE_RV_Results <- read.csv(paste0("/data/williamsjacr/UKB_WGS_Results/BestPRS/RV_",trait,"Best_Betas.csv"))
  RICE_RV_Results$Method <- "RICE-RV"
  
  WGS_Results_Continuous <- rbind(WGS_Results_Continuous,
                                  rbind(CT_Results, LDPred2_Results, LASSOSUM2_Results, RICE_CV_Results, RICE_RV_Results))
}
WGS_Results_Continuous$Data_Type <- "WGS"

## ------------------------------------------------------------
## Theme + colors
## ------------------------------------------------------------
theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) +
      theme(
        plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold", size = 14),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(colour="black", size=2),
        axis.ticks = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.title = element_text(face="bold.italic", size=18),
        legend.text = element_text(size=8),
        plot.margin = unit(c(10,5,5,5), "mm"),
        strip.background = element_rect(colour="#f0f0f0", fill="#f0f0f0"),
        strip.text = element_text(face="bold")
      ))
}

## Named colors (no WES)
scale_fill_Publication <- function(...){
  ggplot2::scale_fill_manual(
    values = c(
      "RICE-CV & Imputed + WES" = "#fab66a",
      "RICE-CV & WGS"           = "#f7911f",
      "RICE-RV & Imputed + WES" = "#3d80ad",
      "RICE-RV & WGS"           = "#16669d"
    ),
    ...
  )
}
scale_fill_Publication_tmp <- scale_fill_Publication

## ------------------------------------------------------------
## Combine + filter (DROP WES here)
## ------------------------------------------------------------
full_results_Continuous <- rbind(Imputed_Results_Continuous, WGS_Results_Continuous)
full_results_Continuous$Method[full_results_Continuous$Method == "CV"] <- "RICE-CV"
full_results_Continuous$Method[full_results_Continuous$Method == "RV"] <- "RICE-RV"
full_results_Continuous <- full_results_Continuous[full_results_Continuous$Method %in% c("RICE-CV","RICE-RV"),]
full_results_Continuous$Method_DataSource <- paste0(full_results_Continuous$Method," & ",full_results_Continuous$Data_Type)

full_results_Continuous <- full_results_Continuous[full_results_Continuous$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_Continuous$trait[full_results_Continuous$trait == "logTG"] <- "log(TG)"
full_results_Continuous$trait <- factor(full_results_Continuous$trait, levels = c("BMI","Height","HDL","LDL","log(TG)","TC"))

full_results_Binary <- rbind(Imputed_Results_Binary, WGS_Results_Binary)
full_results_Binary$Method[full_results_Binary$Method == "CV"] <- "RICE-CV"
full_results_Binary$Method[full_results_Binary$Method == "RV"] <- "RICE-RV"
full_results_Binary <- full_results_Binary[full_results_Binary$Method %in% c("RICE-CV","RICE-RV"),]
full_results_Binary$Method_DataSource <- paste0(full_results_Binary$Method," & ",full_results_Binary$Data_Type)

full_results_Binary <- full_results_Binary[full_results_Binary$ancestry %in% c("AFR","EUR","SAS","AMR"),]
full_results_Binary$trait <- factor(full_results_Binary$trait, levels = c("Asthma","Breast","CAD","Prostate","T2D"))

## keep only these 4 legend entries
mds_levels_main <- c(
  "RICE-CV & Imputed + WES",
  "RICE-CV & WGS",
  "RICE-RV & Imputed + WES",
  "RICE-RV & WGS"
)
full_results_Continuous$Method_DataSource <- factor(full_results_Continuous$Method_DataSource, levels = mds_levels_main)
full_results_Binary$Method_DataSource <- factor(full_results_Binary$Method_DataSource, levels = mds_levels_main)

## clip negatives
full_results_Continuous$beta_adjusted[full_results_Continuous$beta_adjusted < 0] <- 0
full_results_Binary$beta_adjusted[full_results_Binary$beta_adjusted < 0] <- 0

## legend extraction ordering (alternate within dataset)
full_results_tmp <- full_results_Continuous
full_results_tmp$Method_DataSource <- factor(
  full_results_tmp$Method_DataSource,
  levels = c(
    "RICE-CV & Imputed + WES",
    "RICE-RV & Imputed + WES",
    "RICE-CV & WGS",
    "RICE-RV & WGS"
  )
)

ylim_continuous <- max(full_results_Continuous$beta_adjusted, na.rm = TRUE) + 0.03
ylim_binary <- max(full_results_Binary$beta_adjusted, na.rm = TRUE) + 0.03

## ============================================================
## EUR: Continuous + Binary (with shared legend)
## ============================================================
plot1 <- ggplot(full_results_Continuous[full_results_Continuous$ancestry == "EUR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted), fill=Method_DataSource),
           position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0, ylim_continuous) +
  theme_Publication() +
  scale_fill_Publication() +
  guides(fill = guide_legend(title = "PRS Method & Dataset"))

plot_tmp <- ggplot(full_results_tmp[full_results_tmp$ancestry == "EUR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted), fill=Method_DataSource),
           position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Beta of PRS per SD") +
  ylim(0, ylim_continuous) +
  theme_Publication() +
  scale_fill_Publication_tmp() +
  guides(fill = guide_legend(title = "PRS Method & Dataset"))

plot2 <- ggplot(full_results_Binary[full_results_Binary$ancestry == "EUR",]) +
  geom_bar(aes(x=Method, y=abs(beta_adjusted), fill=Method_DataSource),
           position = "dodge", stat="identity", alpha=0.7) +
  facet_grid(cols = vars(trait)) +
  ylab("Log Odds Ratio of PRS per SD") +
  ylim(0, ylim_binary) +
  theme_Publication() +
  scale_fill_Publication() +
  guides(fill = guide_legend(title = "PRS Method & Dataset"))

prow <- plot_grid(
  NULL,
  plot1 + theme(legend.position="none") +
    ggtitle("Comparison of RICE PRS Results using Imputed + WES and WGS for European Ancestry"),
  NULL,
  plot2 + theme(legend.position="none") + theme(plot.title = element_blank()),
  rel_heights = c(-.05, 1, -0.03, 1),
  ncol = 1
)

grob_tmp <- ggplotGrob(plot_tmp)
legend_b <- grob_tmp$grobs[[which(sapply(grob_tmp$grobs, function(x) x$name) == "guide-box")]]

pdf("UKB_ImputedWES_vs_WGS_EUR.pdf", width=10, height=6.18047)
print(plot_grid(prow, NULL, legend_b, NULL, ncol = 1, rel_heights = c(1, -0.02, .1, 0.01)))
dev.off()

## ============================================================
## AFR / AMR / SAS: Continuous only (with legend)
## ============================================================
make_continuous_pdf <- function(anc, out_pdf, title_txt){
  
  p1 <- ggplot(full_results_Continuous[full_results_Continuous$ancestry == anc,]) +
    geom_bar(aes(x=Method, y=abs(beta_adjusted), fill=Method_DataSource),
             position = "dodge", stat="identity", alpha=0.7) +
    facet_grid(cols = vars(trait)) +
    ylab("Beta of PRS per SD") +
    ylim(0, ylim_continuous) +
    theme_Publication() +
    scale_fill_Publication() +
    guides(fill = guide_legend(title="PRS Method & Dataset"))
  
  pt <- ggplot(full_results_tmp[full_results_tmp$ancestry == anc,]) +
    geom_bar(aes(x=Method, y=abs(beta_adjusted), fill=Method_DataSource),
             position = "dodge", stat="identity", alpha=0.7) +
    facet_grid(cols = vars(trait)) +
    ylab("Beta of PRS per SD") +
    ylim(0, ylim_continuous) +
    theme_Publication() +
    scale_fill_Publication_tmp() +
    guides(fill = guide_legend(title="PRS Method & Dataset"))
  
  prow_local <- plot_grid(
    NULL,
    p1 + theme(legend.position="none") + ggtitle(title_txt),
    rel_heights = c(-.05, 1),
    ncol = 1
  )
  
  grob_pt <- ggplotGrob(pt)
  leg <- grob_pt$grobs[[which(sapply(grob_pt$grobs, function(x) x$name) == "guide-box")]]
  
  pdf(out_pdf, width=10, height=6.18047)
  print(plot_grid(prow_local, NULL, leg, NULL, ncol = 1, rel_heights = c(1, -0.02, .1, 0.01)))
  dev.off()
}

make_continuous_pdf(
  anc = "AFR",
  out_pdf = "UKB_ImputedWES_vs_WGS_AFR.pdf",
  title_txt = "Comparison of RICE PRS Results using Imputed + WES and WGS for African Ancestry"
)

make_continuous_pdf(
  anc = "AMR",
  out_pdf = "UKB_ImputedWES_vs_WGS_AMR.pdf",
  title_txt = "Comparison of RICE PRS Results using Imputed + WES and WGS for Admixed American Ancestry"
)

make_continuous_pdf(
  anc = "SAS",
  out_pdf = "UKB_ImputedWES_vs_WGS_SAS.pdf",
  title_txt = "Comparison of RICE PRS Results using Imputed + WES and WGS for South Asian Ancestry"
)
