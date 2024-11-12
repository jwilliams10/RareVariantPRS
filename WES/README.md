## UKB Whole Exome Sequencing Analysis

This folder contains scripts to both clean the UKB WES genotype data and to perform RICE on 11 traits (6 continuous and 5 binary).

### UKB_QC_FullData

This folder contains scripts to clean the raw vcf.gz files. Cleaned versions of these files are used to create the GDS and plink files used for analysis.

### UKB_gds_processing

This folder uses the cleaned vcf.gz files to create a set of gds files to be used by STAARpipeline for rare variant association testing and PRS.

### Common_Variants.sh/Merge_Bed.R, Raw_Phenotypes.R, and Train_Tune_Validation.R

These files create a set of common variants only plink files/merge them into one large file, clean the raw phenotypes, and create train/tune/validation sets of the phenotype data, respectively.

### Binary

Contains all files to perform RICE on the 5 binary traits for UKB WES. Also contains code to create Supplementary Figures 3, 4, and 7.

### Continuous 

Contains all files to perform RICE on the 6 continuous traits for UKB WES. Also contains code to create Figure 4 and Supplementary Figures 4 and 8.

### Remaining Scripts (SampleSizes.R/.sh and RareVariantEffectSizes.R)

The remaining scripts are used to create Supplementary Table 2 and Supplementary Data tables 1-11.