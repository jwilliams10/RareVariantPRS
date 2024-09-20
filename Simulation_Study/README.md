## Simulation_Study

This folder contains code to obtain the burdens and common variants from chromosome 22 of the UKB WES data, generate the simulated outcomes, generate the train/tune/validation splits, and perform RICE for simulated data with 98,343 individuals of EUR ancestry in the training dataset. 

### Full_Data

Contains two folders, Extract_chr22 and G_star. The first folder, Extract_chr22, contains scripts to perform the extraction of chromosome 22 from the UKB WES data. The second folder, G_star, contains functions to call using source() to perform the extraction.

### Simulate_Data

The script inside this folder uses the extracted chromosome 22 data to generate simulated outcomes.

### Train_Tune_Validation_Split

Splits the simulated data into train/tune/validation sets based on the desired training fraction.

### Remaining Scripts

The remaining scripts are almost a carbon copy of the scripts in WES/Continuous/ but slimmed down for speed. The scripts perform RICE for the largest training fraction and are executed using the all.sh script.