## Binary

This folder contains all files to perform RICE on the 5 binary traits for UKB WES. Also contains code to create Supplementary Figure 3, 4, and 7. 

### REGENIE_Reorganize.R/GWAS_SumStats_REGENIE.sh

REGENIE requires a different format of the phenotype file then plink, so REGENIE_Reorganize reorganizes the phenotypes and GWAS_SumStats_REGENIE computes summary statistics for the common variants.

### CommonVariantPRS

This folder contains scripts to perform CT, LDpred2, Lassosum2, and RICE-CV (OneCommonPRS.R). Further contains a script to build the score file for RICE-CV to submit to PGS catalog.

### RareVariant_Analysis

This folder contains scripts to perform gene-centric coding analysis using the STAARpipeline conditional on RICE-CV. First null models are computed for each trait (Null_Models.R), then gene-centric coding analysis is done for regular masks and long masks (2 scripts for each trait submitted with the .sh files), lastly the results are summarized in Modified_Summary_Script.R.

### RareVariant_PRS

This folder contains a script to perform RICE-RV and to obtain coefficients for the significant gene burdens from RICE-RV.

### Common_Plus_Rare_PRS.R/.sh

Performs RICE with RICE-CV and RICE-RV.

### LDSC.R/.sh

Runs LDSC on the regenie summary statistics.

### Overall_Results.R and QQPlots_CV.R

Creates Supplementary Figure 3-4 and Supplementary Figure 7, respectively.