## Integrating Common and Rare Variants Improves Polygenic Risk Prediction Across Diverse Ancestries 

This is a repository containing all analyses conducted for the manuscript, Integrating Common and Rare Variants Improves Polygenic Risk Prediction Across Diverse Ancestries. The analyses are split into 4 components; UKB whole exome sequencing, UKB whole genome sequencing, simulations, and All of Us.

### UKB Whole Exome Sequencing

Scripts to conduct the analysis of UKB Whole Exome Sequencing are in the WES folder.

### UKB Whole Genome Sequencing

Scripts to conduct the analysis of UKB Whole Genome Sequencing are in the DNANexus folder.

### All of Us 

Scripts to conduct the analysis for the All of Us cohort are located in the AoU folder.

### Simulations

Code to simulate the data and to conduct the simulations were split over two folders, Simulation_Study and Simulation_Study2. Simulation_Study contains code to obtain the burdens and common variants from chromosome 22 of the UKB WES data, generate the simulated outcomes, generate the train/tune/validation splits, and perform RICE for simulated data with 98,343 individuals of EUR ancestry in the training dataset. Simulation_Study2 performs RICE for simulated data with 49,173 individuals of EUR ancestry in the training dataset.

### Remaining Scripts

The remaining scripts are to generate figures and tables for the manuscript. Simulation_Results.R creates figures ... Beta_Figure.R creates ... Beta_Table creates Table ... Fig2.R creates ..., Fig5_WES.R creates ... Fig5_WGS.R creates ... QQPlots.R create ... WES_vs_WGS_Results.R creates ... WGS_Ancestry_Adjusted_Plot.R creates ...