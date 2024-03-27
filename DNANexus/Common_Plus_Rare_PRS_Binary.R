rm(list = ls())
library(caret)
library(ranger)
library(SuperLearner)
library(dplyr)
library(boot)
library(RISCA)

# dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/Common_Plus_Rare_PRS_Binary.R -iin=UKB_PRS:JW/UKB_Phenotypes/Scripts/Binary/Common_Plus_Rare_PRS_Binary.sh -icmd="bash Common_Plus_Rare_PRS_Binary.sh" -y --destination UKB_PRS:JW/UKB_Phenotypes/Results/Binary/BestPRS/ --priority high --instance-type mem1_ssd1_v2_x4

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  
  ## STAARO
  
  ## Pull in Phenotypes/Covariates 
  pheno_tuning <- read.delim("All_Tune.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Tune_All.txt"))
  
  rarevariant_prs <- read.delim(paste0(trait,"_Best_All_STAARO_Tune_All.txt"))
  
  file.remove(paste0(trait,"_Best_All_STAARO_Tune_All.txt"))
  
  pheno_tuning <- left_join(pheno_tuning,common_prs,by = "IID")
  colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:26],"CV_PRS")
  pheno_tuning <- left_join(pheno_tuning,rarevariant_prs,by = "IID")
  colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:27],"RV_PRS")
  
  PRSs_Tune <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(27,28)]
  
  if(trait %in% c("Breast","Prostate")){
    confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    tune_model <- glm(as.formula(paste0(trait,"~CV_PRS+RV_PRS+age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data = pheno_tuning,family = binomial)
  }else{
    confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    tune_model <- glm(as.formula(paste0(trait,"~CV_PRS+RV_PRS+age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data = pheno_tuning,family = binomial)
  }
  
  if(is.na(coef(tune_model)[3])){
    pheno_tuning$PRS <- coef(tune_model)[2]*pheno_tuning$CV_PRS
  }else{
    pheno_tuning$PRS <- coef(tune_model)[2]*pheno_tuning$CV_PRS + coef(tune_model)[3]*pheno_tuning$RV_PRS 
  }
  
  
  roc_obj_comb <- roc.binary(status = trait,
                             variable = "PRS",
                             confounders = confounders,
                             data = pheno_tuning[!is.na(pheno_tuning[,trait]),],
                             precision=seq(0.05,0.95, by=0.05))
  
  roc_obj_CV <- roc.binary(status = trait,
                           variable = "CV_PRS",
                           confounders = confounders,
                           data = pheno_tuning[!is.na(pheno_tuning[,trait]),],
                           precision=seq(0.05,0.95, by=0.05))
  
  roc_obj_RV <- roc.binary(status = trait,
                           variable = "RV_PRS",
                           confounders = confounders,
                           data = pheno_tuning[!is.na(pheno_tuning[,trait]),],
                           precision=seq(0.05,0.95, by=0.05))
  
  var <- c("PRS","CV_PRS","RV_PRS")[which.max(c(roc_obj_comb$auc,roc_obj_CV$auc,roc_obj_RV$auc))]
  
  
  
  
  
  pheno_vad <- read.delim("All_Validation.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Validation_All.txt"))
  
  rarevariant_prs <- read.delim(paste0(trait,"_Best_All_STAARO_Validation_All.txt"))
  
  file.remove(paste0(trait,"_Best_All_STAARO_Validation_All.txt"))
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:26],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:27],"RV_PRS")
  
  pheno_vad <- pheno_vad[!is.na(pheno_vad[,trait]),]
  PRSs_Validation <- pheno_vad[,c(1,27,28)]
  
  if(is.na(coef(tune_model)[3])){
    pheno_vad$PRS <- coef(tune_model)[2]*pheno_vad$CV_PRS
  }else{
    pheno_vad$PRS <- coef(tune_model)[2]*pheno_vad$CV_PRS + coef(tune_model)[3]*pheno_vad$RV_PRS 
  }
  
  load("all_phenotypes.RData")
  
  pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  PRSs_Validation_EUR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  PRSs_Validation_NonEur <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  PRSs_Validation_UNK <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  PRSs_Validation_SAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  PRSs_Validation_MIX <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  PRSs_Validation_AFR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  PRSs_Validation_EAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  #EUR
  d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_STAARO_EUR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_EUR.RData")) 
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_STAARO_NonEur",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_NonEur.RData"))  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_STAARO_UNK",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_UNK.RData"))  
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_STAARO_SAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_SAS.RData"))  
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_STAARO_MIX",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_MIX.RData"))  
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_STAARO_AFR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_AFR.RData")) 
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  if(trait %in% c("Prostate","CAD")){
    SL.result <- NA
  }else{
    d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
    
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d,
                          precision=seq(0.05,0.95, by=0.05))
    AUC <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = var,
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
    ci_result <- boot.ci(boot_AUC, type = "perc")
    if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
    
    SL.result <- data.frame(method = "CV_plus_RV_STAARO_EAS",
                            AUC = AUC,
                            AUC_low = ci_result$percent[4],
                            AUC_high = ci_result$percent[5]
    )
  }
  save(SL.result, file = paste0(trait,"_STAARO_All_Result_EAS.RData")) 
  
  
  
  
  
  
  
  ## Burden
  
  ## Pull in Phenotypes/Covariates 
  pheno_tuning <- read.delim("All_Tune.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Tune_All.txt"))
  
  rarevariant_prs <- read.delim(paste0(trait,"_Best_All_Burden_Tune_All.txt"))
  
  file.remove(paste0(trait,"_Best_Tune_All.txt"))
  file.remove(paste0(trait,"_Best_All_Burden_Tune_All.txt"))
  
  pheno_tuning <- left_join(pheno_tuning,common_prs,by = "IID")
  colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:26],"CV_PRS")
  pheno_tuning <- left_join(pheno_tuning,rarevariant_prs,by = "IID")
  colnames(pheno_tuning) <- c(colnames(pheno_tuning)[1:27],"RV_PRS")
  
  PRSs_Tune <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(27,28)]
  
  if(trait %in% c("Breast","Prostate")){
    confounders <- as.formula(paste0("~age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    tune_model <- glm(as.formula(paste0(trait,"~CV_PRS+RV_PRS+age+age2+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data = pheno_tuning,family = binomial)
  }else{
    confounders <- as.formula(paste0("~age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10"))
    tune_model <- glm(as.formula(paste0(trait,"~CV_PRS+RV_PRS+age+age2+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10")),data = pheno_tuning,family = binomial)
  }
  
  if(is.na(coef(tune_model)[3])){
    pheno_tuning$PRS <- coef(tune_model)[2]*pheno_tuning$CV_PRS
  }else{
    pheno_tuning$PRS <- coef(tune_model)[2]*pheno_tuning$CV_PRS + coef(tune_model)[3]*pheno_tuning$RV_PRS 
  }
  
  roc_obj_comb <- roc.binary(status = trait,
                             variable = "PRS",
                             confounders = confounders,
                             data = pheno_tuning[!is.na(pheno_tuning[,trait]),],
                             precision=seq(0.05,0.95, by=0.05))
  
  roc_obj_CV <- roc.binary(status = trait,
                           variable = "CV_PRS",
                           confounders = confounders,
                           data = pheno_tuning[!is.na(pheno_tuning[,trait]),],
                           precision=seq(0.05,0.95, by=0.05))
  
  roc_obj_RV <- roc.binary(status = trait,
                           variable = "RV_PRS",
                           confounders = confounders,
                           data = pheno_tuning[!is.na(pheno_tuning[,trait]),],
                           precision=seq(0.05,0.95, by=0.05))
  
  var <- c("PRS","CV_PRS","RV_PRS")[which.max(c(roc_obj_comb$auc,roc_obj_CV$auc,roc_obj_RV$auc))]
  
  
  
  pheno_vad <- read.delim("All_Validation.txt")
  
  common_prs <- read.delim(paste0(trait,"_Best_Validation_All.txt"))
  
  rarevariant_prs <- read.delim(paste0(trait,"_Best_All_Burden_Validation_All.txt"))
  
  file.remove(paste0(trait,"_Best_Validation_All.txt"))
  file.remove(paste0(trait,"_Best_All_Burden_Validation_All.txt"))
  
  pheno_vad <- left_join(pheno_vad,common_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:26],"CV_PRS")
  pheno_vad <- left_join(pheno_vad,rarevariant_prs,by = "IID")
  colnames(pheno_vad) <- c(colnames(pheno_vad)[1:27],"RV_PRS")
  
  pheno_vad <- pheno_vad[!is.na(pheno_vad[,trait]),]
  PRSs_Validation <- pheno_vad[,c(1,27,28)]
  
  if(is.na(coef(tune_model)[3])){
    pheno_vad$PRS <- coef(tune_model)[2]*pheno_vad$CV_PRS
  }else{
    pheno_vad$PRS <- coef(tune_model)[2]*pheno_vad$CV_PRS + coef(tune_model)[3]*pheno_vad$RV_PRS 
  }
  
  
  load("all_phenotypes.RData")
  
  pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  PRSs_Validation_EUR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  PRSs_Validation_NonEur <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  PRSs_Validation_UNK <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  PRSs_Validation_SAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  PRSs_Validation_MIX <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  PRSs_Validation_AFR <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  PRSs_Validation_EAS <- PRSs_Validation[PRSs_Validation$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  #EUR
  d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_Burden_EUR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_Burden_All_Result_EUR.RData")) 
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_Burden_NonEur",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_Burden_All_Result_NonEur.RData"))  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_Burden_UNK",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_Burden_All_Result_UNK.RData"))  
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_Burden_SAS",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_Burden_All_Result_SAS.RData"))  
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_Burden_MIX",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_Burden_All_Result_MIX.RData"))  
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
  
  roc_obj <- roc.binary(status = trait,
                        variable = var,
                        confounders = confounders,
                        data = d,
                        precision=seq(0.05,0.95, by=0.05))
  AUC <- roc_obj$auc
  
  calc_auc <- function(data, indices) {
    d_sub <- data[indices,] # allows boot to select sample
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d_sub,
                          precision=seq(0.05,0.95, by=0.05))
    return(roc_obj$auc)
  }
  boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
  ci_result <- boot.ci(boot_AUC, type = "perc")
  if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
  
  SL.result <- data.frame(method = "CV_plus_RV_Burden_AFR",
                          AUC = AUC,
                          AUC_low = ci_result$percent[4],
                          AUC_high = ci_result$percent[5]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  save(SL.result, file = paste0(trait,"_Burden_All_Result_AFR.RData")) 
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  if(trait %in% c("Prostate","CAD")){
    SL.result <- NA
  }else{
    d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","PRS","CV_PRS","RV_PRS")]
    
    roc_obj <- roc.binary(status = trait,
                          variable = var,
                          confounders = confounders,
                          data = d,
                          precision=seq(0.05,0.95, by=0.05))
    AUC <- roc_obj$auc
    
    calc_auc <- function(data, indices) {
      d_sub <- data[indices,] # allows boot to select sample
      roc_obj <- roc.binary(status = trait,
                            variable = var,
                            confounders = confounders,
                            data = d_sub,
                            precision=seq(0.05,0.95, by=0.05))
      return(roc_obj$auc)
    }
    boot_AUC <- boot(data = d, statistic = calc_auc, R = 1000)
    ci_result <- boot.ci(boot_AUC, type = "perc")
    if(is.null(ci_result)){ci_result <- data.frame(percent = c(0,0,0,.5,.5))}
    
    SL.result <- data.frame(method = "CV_plus_RV_Burden_EAS",
                            AUC = AUC,
                            AUC_low = ci_result$percent[4],
                            AUC_high = ci_result$percent[5]
    )
  }
  save(SL.result, file = paste0(trait,"_Burden_All_Result_EAS.RData")) 
  
}

system("rm All_Tune.txt")
system("rm All_Validation.txt")
system("rm all_phenotypes.RData")