## Continuous

This folder contains all files to perform RICE on the 6 continuous traits for UKB WES. Also contains code to create Figure 4 and Supplementary Figures 4 and 8.

### GWAS_SumStats.R/.sh

Runs GWAS analysis with plink for the 6 traits.

### CommonVariantPRS

This folder contains scripts to perform CT, LDpred2, Lassosum2, and RICE-CV (OneCommonPRS.R). Further contains a script to build the score file for RICE-CV to submit to PGS catalog.

### RareVariant_Analysis

This folder contains scripts to perform gene-centric coding analysis using the STAARpipeline conditional on RICE-CV. First null models are computed for each trait (Null_Models.R), then gene-centric coding analysis is done for regular masks and long masks (2 scripts for each trait submitted with the .sh files), lastly the results are summarized in Modified_Summary_Script.R.

### RareVariant_PRS

This folder contains a script to perform RICE-RV and to obtain coefficients for the significant gene burdens from RICE-RV.

### Common_Plus_Rare_PRS.R/.sh

Performs RICE with RICE-CV and RICE-RV.

### Overall_Results.R and QQPlots_CV.R

Creates Figure 4/Supplementary Figure 4 and Supplementary Figure 8, respectively.