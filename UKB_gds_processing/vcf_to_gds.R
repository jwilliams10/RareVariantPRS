print(commandArgs(TRUE))
chr <- as.numeric(commandArgs(TRUE)[1])
print(chr)

library(SeqArray)

if(!file.exists(paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/all_chr",chr,".gds"))){
  SeqArray::seqVCF2GDS(vcf.fn = paste0("/data/williamsjacr/UKB_WES_lipids/Data/pVCF/chr",chr,"/ukbb_wes_200k_chr",chr,".vcf.bgz"),out.fn = paste0("/data/williamsjacr/UKB_WES_lipids/Data/gds/all_chr",chr,".gds"),storage.option="LZMA_RA")  
}