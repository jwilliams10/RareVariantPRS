# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Train_Tune_Validation.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Train_Tune_Validation.sh  -icmd="bash Train_Tune_Validation.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/ --instance-type mem1_ssd1_v2_x4

print("ls:")
system("ls")

if(!("data.table" %in% rownames(installed.packages()))){
  install.packages("data.table",quiet = TRUE)
}

if(!("dplyr" %in% rownames(installed.packages()))){
  install.packages("dplyr",quiet = TRUE)
}

if(!("SeqArray" %in% rownames(installed.packages()))){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",quiet = TRUE)
  
  BiocManager::install("SeqArray")
}

library(dplyr)
library(data.table)
library(SeqArray)

ukb_pheno <- fread("ukb_multi.pheno")
ukb_covariate <- fread("ukb_multi.cov")
ukb_ancestries <- readRDS("ukb_multi_anc.RDS")




ukb_covariate <- ukb_covariate[,c("FID","IID","ancestry","female","age",paste0("pc",1:10))]
colnames(ukb_covariate) <- c(colnames(ukb_covariate)[1:3],"sex",colnames(ukb_covariate)[5:15])

ukb_pheno <- inner_join(ukb_pheno,ukb_covariate)

colnames(ukb_pheno)[colnames(ukb_pheno) == "ancestry"] <- "ethnicity"

ukb_pheno <- inner_join(ukb_pheno,ukb_ancestries[,c("FID","predicted")])
colnames(ukb_pheno)[colnames(ukb_pheno) == "predicted"] <- "ancestry"

ukb_pheno <- as.data.frame(ukb_pheno)

bed_ids <- fread("all_chr.fam")
bed_ids <- as.data.frame(bed_ids)
bed_ids <- bed_ids[,1:2]

gds <- seqOpen("ukb.200k.wgs.chr22.pass.annotated.gds")
gds_ids <- seqGetData(gds, "sample.id")
seqClose(gds)

sampleids <- unique(c(bed_ids[,1],gds_ids))

sampleids <- bed_ids[bed_ids[,1] %in% sampleids,]

ukb_pheno <- ukb_pheno[ukb_pheno$IID %in% sampleids[,1],]
sampleids <- sampleids[sampleids[,1] %in% ukb_pheno$IID,]

sampleids <- sampleids[order(sampleids[,1]),,drop = FALSE]
ukb_pheno <- ukb_pheno[order(ukb_pheno$IID),]

set.seed(1335)

i <- (1:nrow(sampleids))[ukb_pheno$ancestry == "EUR"]

train_number <- round(nrow(sampleids)*0.7) + 1
train <- sample(i, train_number)

i <- (1:nrow(sampleids))[!((1:nrow(sampleids)) %in% train)]
i_EUR <- i[ukb_pheno$ancestry[i] == "EUR"]
i_AFR <- i[ukb_pheno$ancestry[i] == "AFR"]
i_SAS <- i[ukb_pheno$ancestry[i] == "SAS"]
i_EAS <- i[ukb_pheno$ancestry[i] == "EAS"]
i_AMR <- i[ukb_pheno$ancestry[i] == "AMR"]

tune <- c(sample(i_EUR,round(length(i_EUR)/2)),
          sample(i_AFR,round(length(i_AFR)/2)),
          sample(i_SAS,round(length(i_SAS)/2)),
          sample(i_EAS,round(length(i_EAS)/2)),
          sample(i_AMR,round(length(i_AMR)/2)))

validation <- i[!(i %in% tune)]

train <- sampleids[train,]
tune <- sampleids[tune,]
validation <- sampleids[validation,]

write.table(train,"train.txt",row.names = FALSE,col.names = FALSE)
write.table(tune,"tune.txt",row.names = FALSE,col.names = FALSE)
write.table(validation,"validation.txt",row.names = FALSE,col.names = FALSE)

reference <- train[sample(1:nrow(train),3000,replace = FALSE),]

write.table(reference,"reference.txt",row.names = FALSE,col.names = FALSE)

save(ukb_pheno,file = "all_phenotypes.RData")


nrow(train)
nrow(tune)
nrow(validation)


phenotype_train <- ukb_pheno[ukb_pheno$IID %in% train[,1],]
phenotype_tune <- ukb_pheno[ukb_pheno$IID %in% tune[,1],]
phenotype_validation <- ukb_pheno[ukb_pheno$IID %in% validation[,1],]

# phenotype_train$FID <- 0
# phenotype_tune$FID <- 0
# phenotype_validation$FID <- 0

phenotype_train <- phenotype_train[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
phenotype_train$age2 <- phenotype_train$age^2
phenotype_tune <- phenotype_tune[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
phenotype_tune$age2 <- phenotype_tune$age^2
phenotype_validation <- phenotype_validation[,c("IID","FID","BMI","TC","HDL","LDL","logTG","Height","Asthma","CAD","T2D","Breast","Prostate","age","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")]
phenotype_validation$age2 <- phenotype_validation$age^2

colnames(phenotype_train) <- c("IID",colnames(phenotype_train)[-1])
colnames(phenotype_tune) <- c("IID",colnames(phenotype_tune)[-1])
colnames(phenotype_validation) <- c("IID",colnames(phenotype_validation)[-1])

phenotype_train[c(14,16:26)] <- lapply(phenotype_train[c(14,16:26)], function(x){c(scale(x))})
phenotype_tune[c(14,16:26)] <- lapply(phenotype_tune[c(14,16:26)], function(x){c(scale(x))})
phenotype_validation[c(14,16:26)] <- lapply(phenotype_validation[c(14,16:26)], function(x){c(scale(x))})

save(phenotype_train,file = "All_Train.RData")
save(phenotype_tune,file = "All_Tune.RData")
save(phenotype_validation,file = "All_Validation.RData")

write.table(phenotype_train,file = "All_Train.txt",sep = '\t',row.names = FALSE,quote = FALSE)
write.table(phenotype_tune,file = "All_Tune.txt",sep = '\t',row.names = FALSE,quote = FALSE)
write.table(phenotype_validation,file = "All_Validation.txt",sep = '\t',row.names = FALSE,quote = FALSE)

file.remove("all_chr.fam")
file.remove("ukb.200k.wgs.chr22.pass.annotated.gds")
file.remove("ukb_multi.cov")
file.remove("ukb_multi.pheno")
file.remove("ukb_multi_anc.RDS")
