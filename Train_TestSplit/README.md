SampleIds.sh creates a file of all sample ids in chr1.
Train_Tune_Validate_Split.R splits the sample ids to train, test, validate sets and stores the sample ids for each of those sets after filtering for LDL missingness and unrelatedness.
The sample ids created in SampleIDs.sh contains 14 negative id values. When obtaining the subsets these values are then lost in which ever subset they appeared in.
RareVariants_TrainTuneValidation.sh takes the full data (common and rare) and builds the respective subsets as vcf.bgz files.
CommonVariants_TrainTuneValidation.sh first creates a full common variant .pgen file using MAF 0.01 filter with plink2. The creates the respective subsets as .pgen files.