## Continuous

This folder contains all files to perform RICE on the 6 continuous traits for UKB WES. Also contains code to create Supplementary Figures 3, 4, and 7.

### GWAS_SumStats_REGENIE.sh

Runs GWAS analysis with regenie for the 6 traits.

### CommonVariantPRS

This folder contains scripts to perform CT, LDpred2, Lassosum2, and RICE-CV (OneCommonPRS.R). Also contains a script to build the score file for RICE-CV to submit to PGS catalog.

### RareVariant_Analysis

This folder contains scripts to perform gene-centric coding analysis using the STAARpipeline conditional on RICE-CV. First null models are computed for each trait (Null_Models.R), then gene-centric coding analysis is done for regular masks and long masks (2 scripts for each trait submitted with the .sh files), lastly the results are summarized in Modified_Summary_Script.R.

### RareVariant_PRS

This folder contains a script to perform RICE-RV and to obtain coefficients for the significant gene burdens from RICE-RV.

### Common_Plus_Rare_PRS.R/.sh

Performs RICE with RICE-CV and RICE-RV.

### Overall_Results.R and QQPlots_CV_REGENIE.R

Creates Supplementary Figures 3/4 and Supplementary Figure 7, respectively.

### LDSC.R/.sh

Runs LDSC on the regenie summary statistics.

### All_Train_LDSC.R, GWAS_SumStats_LDSC.R/.sh, and LDSC_Results.R

Sensitivity analysis using different GWAS methods and responses.