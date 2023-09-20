print(commandArgs(TRUE))
chr <- as.numeric(commandArgs(TRUE)[1])

library(SeqArray)

if(!file.exists(paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/train_chr",chr,".gds"))){
  SeqArray::seqVCF2GDS(paste0("/data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr",chr,"/ukbb_wes_200k_chr",chr,"_train.vcf.bgz"), out.fn = paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/train_chr",chr,".gds"), header = NULL, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)
}
if(!file.exists(paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/tune_chr",chr,".gds"))){
  SeqArray::seqVCF2GDS(paste0("/data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr",chr,"/ukbb_wes_200k_chr",chr,"_tune.vcf.bgz"), out.fn = paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/tune_chr",chr,".gds"), header = NULL, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)
}
if(!file.exists(paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/validation_chr",chr,".gds"))){
  SeqArray::seqVCF2GDS(paste0("/data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr",chr,"/ukbb_wes_200k_chr",chr,"_validation.vcf.bgz"), out.fn = paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/validation_chr",chr,".gds"), header = NULL, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)
}