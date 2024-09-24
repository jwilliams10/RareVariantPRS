rm(list = ls())

load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Train.RData")
load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Tune.RData")
load("/data/williamsjacr/UKB_WES_Simulation/Simulation1/simulated_data/phenotypes/Y_Validation.RData")

pheno_train <- Y_train[[1]]
pheno_tuning <- Y_tune[[1]]
pheno_vad <- Y_validation[[1]]

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_train_EUR <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_train_SAS <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_train_AMR <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_train_AFR <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_train_EAS <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_tuning_EUR <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_tuning_SAS <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_tuning_AMR <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_tuning_AFR <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_tuning_EAS <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_vad_EUR <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_AMR <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

samps <- data.frame(Training_Fraction = nrow(pheno_train),Ancestry = c("EUR","SAS","AMR","AFR","EAS"),
                    Training_Sample_Size = c(nrow(pheno_train_EUR),nrow(pheno_train_SAS),nrow(pheno_train_AMR),nrow(pheno_train_AFR),nrow(pheno_train_EAS)),
                    Tuning_Sample_Size = c(nrow(pheno_tuning_EUR),nrow(pheno_tuning_SAS),nrow(pheno_tuning_AMR),nrow(pheno_tuning_AFR),nrow(pheno_tuning_EAS)),
                    Validation_Sample_Size = c(nrow(pheno_vad_EUR),nrow(pheno_vad_SAS),nrow(pheno_vad_AMR),nrow(pheno_vad_AFR),nrow(pheno_vad_EAS)))

load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Train.RData")
load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Tune.RData")
load("/data/williamsjacr/UKB_WES_Simulation/Simulation2/simulated_data/phenotypes/Y_Validation.RData")

pheno_train <- Y_train[[1]]
pheno_tuning <- Y_tune[[1]]
pheno_vad <- Y_validation[[1]]

load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")

pheno_train_EUR <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_train_SAS <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_train_AMR <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_train_AFR <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_train_EAS <- pheno_train[pheno_train$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_tuning_EUR <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_tuning_SAS <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_tuning_AMR <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_tuning_AFR <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_tuning_EAS <- pheno_tuning[pheno_tuning$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

pheno_vad_EUR <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
pheno_vad_SAS <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
pheno_vad_AMR <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AMR"],]
pheno_vad_AFR <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
pheno_vad_EAS <- pheno_vad[pheno_vad$IDs %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]

samps <- rbind(samps,data.frame(Training_Fraction = nrow(pheno_train),Ancestry = c("EUR","SAS","AMR","AFR","EAS"),
                    Training_Sample_Size = c(nrow(pheno_train_EUR),nrow(pheno_train_SAS),nrow(pheno_train_AMR),nrow(pheno_train_AFR),nrow(pheno_train_EAS)),
                    Tuning_Sample_Size = c(nrow(pheno_tuning_EUR),nrow(pheno_tuning_SAS),nrow(pheno_tuning_AMR),nrow(pheno_tuning_AFR),nrow(pheno_tuning_EAS)),
                    Validation_Sample_Size = c(nrow(pheno_vad_EUR),nrow(pheno_vad_SAS),nrow(pheno_vad_AMR),nrow(pheno_vad_AFR),nrow(pheno_vad_EAS))))

samps <- samps[samps$Ancestry != "EAS",]

samps$Training_Fraction <- as.character(samps$Training_Fraction)
samps$Training_Fraction[samps$Training_Fraction == "98343"] <- "98,343"
samps$Training_Fraction[samps$Training_Fraction == "49173"] <- "49,173"
samps$Training_Fraction <- factor(samps$Training_Fraction,levels = c("49,173","98,343"))
samps$Ancestry <- factor(samps$Ancestry,levels = c("AFR","AMR","EUR","SAS"))
samps <- samps[order(samps$Training_Fraction,samps$Ancestry),]
