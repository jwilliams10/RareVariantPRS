## Imputed

This folder contains all files to perform RICE on the UKB imputed genotype + whole exome sequencing as well as additional sensitivity analysis.

### GWAS_SumStats_Continuous.sh and GWAS_SumStats_Binary.sh

Runs GWAS analysis with regenie for the 6 continuous traits and 5 binary traits respectively.

### CommonVariantPRS

This folder contains scripts to perform CT, LDpred2, Lassosum2, and RICE-CV (OneCommonPRS_All.R). Also contains a script to build the score file for RICE-CV to submit to PGS catalog. CT_NewPCs is to a sensitivity analysis.

### RareVariant_Analysis

This folder contains scripts to perform gene-centric coding analysis using the STAARpipeline conditional on RICE-CV. First null models are computed for each trait (Null_Models.R), then gene-centric coding analysis is done for regular masks and long masks, lastly the results are summarized in Modified_Summary_Script.R.

### RareVariant_PRS

This folder contains a script (Single_RareVariant_PRS_All.R/.sh) to perform RICE-RV and to obtain coefficients for the significant gene burdens from RICE-RV. Remaining scripts are to perform sensitivity analyses.

### Common_Plus_Rare_PRS.R/.sh

Performs RICE with RICE-CV and RICE-RV.

### Overall_Results.R and QQPlots_CV.R

Creates Supplementary Figures 9/10 and Supplementary Figure 13, respectively.

### LDSC.R/.sh

Runs LDSC on the regenie summary statistics.

### Scripts with NewPCs in title.

Sensitivity analysis using different PCs.

### Time_Table.R and SampleSizes.R

Used to create Supplementary Table 3 and 12.