## Integrating Common and Rare Variants Improves Polygenic Risk Prediction Across Diverse Ancestries 

This is a repository containing all analyses conducted for the manuscript, Integrating Common and Rare Variants Improves Polygenic Risk Prediction Across Diverse Ancestries. The analyses are split into 5 components; UKB whole exome sequencing, UKB imputed genotype + whole exome sequencing, UKB whole genome sequencing, simulations, and All of Us.

### UKB Whole Exome Sequencing

Scripts to conduct the analysis of UKB Whole Exome Sequencing are in the WES folder.

### UKB imputed genotype + whole exome sequencing

Scripts to conduct the analysis of UKB imputed genotype + whole exome sequencing are in the Imputed folder.

### UKB Whole Genome Sequencing

Scripts to conduct the analysis of UKB Whole Genome Sequencing are in the DNANexus folder.

### All of Us 

Scripts to conduct the analysis for the All of Us cohort are located in the AoU folder.

### Simulations

Code to simulate the data and to conduct the simulations were split across four folders, Simulation_Study, Simulation_Study2, Simulation_Study3, and Simulation_Study4. Simulation_Study contains code to obtain the burdens and common variants from chromosome 22 of the UKB WES data, generate the simulated outcomes, generate the train/tune/validation splits, and perform RICE for simulated data with 98,343 individuals of EUR ancestry in the training dataset. Simulation_Study2 performs RICE for simulated data with 49,173 individuals of EUR ancestry in the training dataset Simulation 3 and 4 are with 98,343 individuals and 49,173 individuals of EUR ancestry in the training dataset, respectively. Simulation 1 and 2 are assuming all rare variants within a causal rare variant set are causal while Simulation 3 and 4 assume a proportion of rare variants within a causal rare variant set are causal.

### Remaining Scripts

The remaining scripts are to generate figures and tables for the manuscript. Simulation_Results.R creates Figure 3 and Supplementary Figure 1; Simulation_SampleSize.R creates Supplementary Table 1; Sim_Characteristics.R creates Supplementary Figure 2; Metrics_Table.R creates Supplementary Table 7, 8, 9, and 10; Fig2.R creates Figure 2; Fig5_WES.R creates Supplementary Figures 5 and 6; Fig5_Imputed.R creates Supplementary Figures 11 and 12; Fig5_WGS.R creates Figure 5 and Supplementary Figures 17 and 18; QQPlots.R creates Supplementary Figures 8, 14, 22, and 27; WES_vs_WGS_Results.R creates Figure 6 and Supplementary Figure 19; RareVariantEffectSizes.R creates the Supplementary Data table; LDSC_Results.R creates Supplementary Table 13.