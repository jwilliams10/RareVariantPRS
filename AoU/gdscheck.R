rm(list = ls())

Sys.getenv("WORKSPACE_BUCKET")
Sys.getenv("OWNER_EMAIL")


BiocManager::install("SeqVarTools")
BiocManager::install("GENESIS")
devtools::install_github("xihaoli/STAAR")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
devtools::install_github("zilinli1988/SCANG")
devtools::install_github("xihaoli/MultiSTAAR",ref="main")
devtools::install_github("xihaoli/STAARpipeline",ref="main")
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Matrix)
library(SCANG)
library(STAARpipeline)



ncRNA <- function (chr, gene_name, genofile, obj_nullmodel, rare_maf_cutoff = 0.01, 
                   rv_num_cutoff = 2, QC_label = "annotation/filter", variant_type = c("SNV","Indel", "variant"), geno_missing_imputation = c("mean","minor"), Annotation_dir = "annotation/info/FunctionalAnnotation", 
                   Annotation_name_catalog, Use_annotation_weights = c(TRUE,FALSE), Annotation_name = NULL, SPA_p_filter = TRUE, 
                   p_filter_cutoff = 0.05, silent = FALSE){
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)
  phenotype.id <- obj_nullmodel$id_include
  n_pheno <- obj_nullmodel$n.pheno
  if (!is.null(obj_nullmodel$use_SPA)) {
    use_SPA <- obj_nullmodel$use_SPA
  }else {
    use_SPA <- FALSE
  }
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- filter == "PASS"
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter)
  gc()
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir, 
                                                  Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                      "GENCODE.Category")]))
  is.in <- ((GENCODE.Category == "ncRNA_exonic") | (GENCODE.Category == 
                                                      "ncRNA_exonic;splicing") | (GENCODE.Category == "ncRNA_splicing")) & 
    (SNVlist)
  variant.id.ncRNA <- variant.id[is.in]
  rm(GENCODE.Category)
  gc()
  seqSetFilter(genofile, variant.id = variant.id.ncRNA, sample.id = phenotype.id)
  rm(variant.id.ncRNA)
  gc()
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir, 
                                              Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                  "GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
  Gene <- as.character(sapply(GENCODE.Info.split, function(z) gsub("\\(.*\\)", 
                                                                   "", z[1])))
  Gene_list_1 <- as.character(sapply(strsplit(Gene, ","), "[", 
                                     1))
  Gene_list_2 <- as.character(sapply(strsplit(Gene, ","), "[", 
                                     2))
  Gene_list_3 <- as.character(sapply(strsplit(Gene, ","), "[", 
                                     3))
  rm(GENCODE.Info)
  gc()
  rm(GENCODE.Info.split)
  gc()
  variant.id.ncRNA <- seqGetData(genofile, "variant.id")
  seqResetFilter(genofile)
  is.in <- union(which(Gene_list_1 == gene_name), which(Gene_list_2 == 
                                                          gene_name))
  is.in <- union(is.in, which(Gene_list_3 == gene_name))
  variant.is.in <- variant.id.ncRNA[is.in]
  seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
  id.genotype <- seqGetData(genofile, "sample.id")
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, 
                                                           length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, 
                                         id.genotype.merge, by = c(phenotype.id = "id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop = FALSE]
  if (!is.null(dim(Geno))) {
    if (dim(Geno)[2] > 0) {
      if (geno_missing_imputation == "mean") {
        Geno <- matrix_flip_mean(Geno)$Geno
      }
      if (geno_missing_imputation == "minor") {
        Geno <- matrix_flip_minor(Geno)$Geno
      }
    }
  }
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                       Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, 
                                                          Annotation_name_catalog$dir[which(Annotation_name_catalog$name == 
                                                                                              Annotation_name[k])]))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, 
                                      Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, 
                                         paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, 
                                      Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  pvalues <- 0
  if (n_pheno == 1) {
    if (!use_SPA) {
      try(pvalues <- STAAR(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                           rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff), 
          silent = silent)
    } else {
      try(pvalues <- STAAR_Binary_SPA(Geno, obj_nullmodel, 
                                      Anno.Int.PHRED.sub, rare_maf_cutoff = rare_maf_cutoff, 
                                      rv_num_cutoff = rv_num_cutoff, SPA_p_filter = SPA_p_filter, 
                                      p_filter_cutoff = p_filter_cutoff), silent = silent)
    }
  } else {
    try(pvalues <- MultiSTAAR(Geno, obj_nullmodel, Anno.Int.PHRED.sub, 
                              rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff), 
        silent = silent)
  }
  results <- c()
  if (inherits(pvalues, "list")) {
    results_temp <- rep(NA, 4)
    results_temp[3] <- "ncRNA"
    results_temp[2] <- chr
    results_temp[1] <- as.character(gene_name)
    results_temp[4] <- pvalues$num_variant
    if (!use_SPA) {
      results_temp <- c(results_temp, pvalues$cMAC, pvalues$results_STAAR_S_1_25, 
                        pvalues$results_STAAR_S_1_1, pvalues$results_STAAR_B_1_25, 
                        pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_A_1_25, 
                        pvalues$results_STAAR_A_1_1, pvalues$results_ACAT_O, 
                        pvalues$results_STAAR_O)
    } else {
      results_temp <- c(results_temp, pvalues$cMAC, pvalues$results_STAAR_B_1_25, 
                        pvalues$results_STAAR_B_1_1, pvalues$results_STAAR_B)
    }
    results <- rbind(results, results_temp)
  }
  if (!is.null(results)) {
    if (!use_SPA) {
      colnames(results) <- colnames(results, do.NULL = FALSE, 
                                    prefix = "col")
      colnames(results)[1:5] <- c("Gene name", "Chr", "Category", 
                                  "#SNV", "cMAC")
      colnames(results)[(dim(results)[2] - 1):dim(results)[2]] <- c("ACAT-O", 
                                                                    "STAAR-O")
    } else {
      colnames(results) <- colnames(results, do.NULL = FALSE, 
                                    prefix = "col")
      colnames(results)[1:5] <- c("Gene name", "Chr", "Category", 
                                  "#SNV", "cMAC")
      colnames(results)[dim(results)[2]] <- c("STAAR-B")
    }
  }
  seqResetFilter(genofile)
  return(results)
}





system("mkdir aouprs")
setwd("/home/rstudio/aouprs/")
system("mkdir dataAGDS")

path <- "gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/dataAGDS"
system(paste0("gsutil ls ",path))

system(paste0("gsutil cp ", path,"/acaf_threshold.chr",21,".gds", " dataAGDS/"))# copy the gds file from google cloud to local workspace

genofile <- seqOpen("dataAGDS/acaf_threshold.chr21.gds")
sampids <- seqGetData(genofile,"sample.id")
varids <- seqGetData(genofile,"variant.id")
varInfo <- seqGetData(genofile,"annotation/info/FunctionalAnnotation/VarInfo")
genofile
seqClose(genofile)


system("gsutil ls gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/Tony_Files")
system("mkdir phenotypes")
system("gsutil cp gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/Tony_Files/aou_BMI_cov.tsv phenotypes/")
system("gsutil cp gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/Tony_Files/aou_BMI_pheno.tsv phenotypes/")


aou_BMI_cov <- read.delim("~/aouprs/phenotypes/aou_BMI_cov.tsv")
aou_BMI_pheno <- read.delim("~/aouprs/phenotypes/aou_BMI_pheno.tsv")
aou_BMI <- inner_join(aou_BMI_pheno,aou_BMI_cov)
aou_BMI$age2 <- aou_BMI$age^2
aou_BMI[c(4,6:26)] <- lapply(aou_BMI[c(4,6:26)], function(x){c(scale(x))})

obj.STAAR.UKB <- fit_nullmodel(as.formula(paste0("BMI","~age+age2+female+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")), data = aou_BMI,id = "IID",kins = NULL,family = gaussian(link = "identity"))

genofile <- seqOpen("dataAGDS/acaf_threshold.chr21.gds", readonly = FALSE)
QC_label <- rep("PASS",length(varids)) 

Anno.folder <- index.gdsn(genofile, "annotation/info")
add.gdsn(Anno.folder, "QC_label", val=QC_label, compress="LZMA_ra", closezip=TRUE)

seqClose(genofile)














## Null model
obj_nullmodel <- obj.STAAR.UKB

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- data.frame(name = c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer",
                                               "CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed",
                                               "aPC.EpigeneticTranscription","aPC.Conservation","aPC.LocalDiversity","aPC.Mappability",
                                               "aPC.TF","aPC.Protein"),
                                      dir = c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info","/genecode_comprehensive_exonic_category",
                                              "/metasvm_pred","/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
                                              "/apc_epigenetics_active","/apc_epigenetics_repressed",
                                              "/apc_epigenetics_transcription","/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
                                              "/apc_transcription_factor","/apc_protein_function"))
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")





## have to make sure this has line 101 being 21
arrayid <- 215

## gene number in job
gene_num_in_array <- 100 
group.num.allchr <- ceiling(table(ncRNA_gene[,1])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
group.num <- group.num.allchr[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
sub_seq_num <- dim(ncRNA_gene_chr)[1]

if(groupid < group.num)
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}

## exclude large ncRNA masks
if(arrayid==117)
{
  sub_seq_id <- setdiff(sub_seq_id,53)
}

if(arrayid==218)
{
  sub_seq_id <- setdiff(sub_seq_id,19)
}

if(arrayid==220)
{
  sub_seq_id <- setdiff(sub_seq_id,c(208,274))
}

if(arrayid==221)
{
  sub_seq_id <- setdiff(sub_seq_id,311)
}

if(arrayid==156)
{
  sub_seq_id <- setdiff(sub_seq_id,41)
}

if(arrayid==219)
{
  sub_seq_id <- setdiff(sub_seq_id,103)
}

## aGDS file
genofile <- seqOpen("dataAGDS/acaf_threshold.chr21.gds")

results_ncRNA <- c()
for(kk in sub_seq_id){
  print(kk)
  gene_name <- ncRNA_gene_chr[kk,2]
  results <- c()
  results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
  
  results_ncRNA <- rbind(results_ncRNA,results)
}

seqClose(genofile)
