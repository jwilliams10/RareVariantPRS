## DNAnexus

This folder contains scripts to perform the analysis of 6 continuous and 5 binary traits using UKB whole genome sequencing (WGS) data on the DNAnexus platform. The outline of this folder is a free-for-all but the scripts roughly follow the scripts in WES/Continuous and WES/Binary with _Binary distinguishing if the script is used for binary traits or not. The scripts that are not the same in either WES/Continuous or WES/Binary are discussed below. The analysis were mostly conducted using the swiss-army-knife tool provided by DNAnexus built using a publicly available docker (https://hub.docker.com/repository/docker/willja16/r_with_plink/general). Some analyses such as REGENIE were ran with tools provided by DNAnexus. Supplementary Figures 10, 15, and 16 are created in Overall_Results.R, Overall_Result_Binary.R, QQPlots_CV.R, QQPlots_CV_Binary.R. 

### Extract_Results.R

Used to extract the results to a zip folder for easier download.

### G_extraction_...

These set of files were used to build the burdens for RICE-RV. Originally this was incorporated with the Single_RareVariant_PRS.R script in WES/Continuous or WES/Binary, but the nature of DNAnexus required this to be done in a separate script.

### Rare Variant Analysis

STAARpipeline is available as both a RAP (https://github.com/li-lab-genetics/staarpipeline-rap) and a R package. To conduct rare variant analysis for gene-centric coding and ncRNA, the RAP was used. For version issues with other packages, gene-centric non-coding analysis was conducted using the R package. The commands to implement the RAP are in the RareVariant_Analysis_Commands text folder and the scripts to implement the gene-centric noncoding analyses is RareVariant_Analysis_Noncoding_Preload.R/.sh.
