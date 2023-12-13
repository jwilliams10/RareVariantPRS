rm(list = ls())
library(data.table)
library(dplyr)
library(RISCA)
library(boot)
library(pROC)

ct_result_rawAUC <- NULL

library(ggplot2)

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_blank(), 
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}  

trait <- "Asthma"

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
  
  prs_mat_train <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_train.txt"), sep="")
  prs_mat_tune <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_tune.txt"), sep="")
  prs_mat_validation <- read.csv(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_all_validation.txt"), sep="")
  
  
  ## Pull in Phenotypes/Covariates 
  pheno_train <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Train.txt")
  pheno_train <- left_join(pheno_train,prs_mat_train,by = "IID")
  
  pheno_tuning <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Tune.txt")
  pheno_tuning <- left_join(pheno_tuning,prs_mat_tune,by = "IID")
  
  pheno_vad <- read.delim("/data/williamsjacr/UKB_WES_Phenotypes/All_Validation.txt")
  pheno_vad <- left_join(pheno_vad,prs_mat_validation,by = "IID")
  
  load("/data/williamsjacr/UKB_WES_Phenotypes/all_phenotypes.RData")
  
  pheno_vad_EUR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EUR"],]
  pheno_vad_NonEur <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry != "EUR"],]
  pheno_vad_UNK <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "UNK"],]
  pheno_vad_SAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "SAS"],]
  pheno_vad_MIX <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "MIX"],]
  pheno_vad_AFR <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "AFR"],]
  pheno_vad_EAS <- pheno_vad[pheno_vad$IID %in% ukb_pheno$IID[ukb_pheno$ancestry == "EAS"],]
  
  #calculate AUC for each of the tuning dataset
  # This is done by regressing the residuals of the model with all covariates against the prs
  AUC_tun_vec <- rep(0,length(pthres))
  
  k <- 1
  for(k in 1:length(pthres)){
    d <- pheno_tuning[!is.na(pheno_tuning[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",k))]
    
    a <- as.formula(paste0(trait,"~p_value_",k))
    roc_obj <- roc(a,data = d,quiet = TRUE)
    AUC_tun_vec[k] <- auc(roc_obj)
  }
  
  #find best p-value threshold
  idx <- which.max(AUC_tun_vec)
  #write down best prs in the prs folder/ save it and store it
  # prs_train_max <- pheno_train[,c("IID","FID",paste0("p_value_",idx))]
  # colnames(prs_train_max) <- c("IID","FID","prs")
  # write.table(prs_train_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_train_best.txt"),row.names = F)
  # 
  # prs_tune_max <- pheno_tuning[,c("IID","FID",paste0("p_value_",idx))]
  # colnames(prs_tune_max) <- c("IID","FID","prs")
  # write.table(prs_tune_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_tune_best.txt"),row.names = F)
  # 
  # prs_vad_max <- pheno_vad[,c("IID","FID",paste0("p_value_",idx))]
  # colnames(prs_vad_max) <- c("IID","FID","prs")
  # write.table(prs_vad_max, file = paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",trait,"_prs_validation_best.txt"),row.names = F)
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_EUR[!is.na(pheno_vad_EUR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  a <- as.formula(paste0(trait,"~p_value_",idx))
  roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
  AUC <- auc(roc_obj)
  ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
  
  AUC.result <- data.frame(method = "CT_EUR",
                           AUC = AUC,
                           AUC_low = ci_result[1],
                           AUC_high = ci_result[3]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result_EUR <- AUC.result
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_NonEur[!is.na(pheno_vad_NonEur[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  a <- as.formula(paste0(trait,"~p_value_",idx))
  roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
  AUC <- auc(roc_obj)
  ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
  
  AUC.result <- data.frame(method = "CT_NonEur",
                           AUC = AUC,
                           AUC_low = ci_result[1],
                           AUC_high = ci_result[3]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result_NonEUR <- AUC.result
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_UNK[!is.na(pheno_vad_UNK[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  a <- as.formula(paste0(trait,"~p_value_",idx))
  roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
  AUC <- auc(roc_obj)
  ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
  
  AUC.result <- data.frame(method = "CT_UNK",
                           AUC = AUC,
                           AUC_low = ci_result[1],
                           AUC_high = ci_result[3]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result_UNK <- AUC.result
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_SAS[!is.na(pheno_vad_SAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  a <- as.formula(paste0(trait,"~p_value_",idx))
  roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
  AUC <- auc(roc_obj)
  ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
  
  AUC.result <- data.frame(method = "CT_SAS",
                           AUC = AUC,
                           AUC_low = ci_result[1],
                           AUC_high = ci_result[3]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result_SAS <- AUC.result
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_MIX[!is.na(pheno_vad_MIX[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  a <- as.formula(paste0(trait,"~p_value_",idx))
  roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
  AUC <- auc(roc_obj)
  ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
  
  AUC.result <- data.frame(method = "CT_MIX",
                           AUC = AUC,
                           AUC_low = ci_result[1],
                           AUC_high = ci_result[3]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result_MIX <- AUC.result
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  d <- pheno_vad_AFR[!is.na(pheno_vad_AFR[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
  
  a <- as.formula(paste0(trait,"~p_value_",idx))
  roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
  AUC <- auc(roc_obj)
  ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
  
  AUC.result <- data.frame(method = "CT_AFR",
                           AUC = AUC,
                           AUC_low = ci_result[1],
                           AUC_high = ci_result[3]
  )
  
  ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
  ct.result_AFR <- AUC.result
  
  
  #evaluate the best threshold based on the tuning on the validation dataset
  if(trait %in% c("Prostate","CAD")){
    AUC.result <- NA
  }else{
    d <- pheno_vad_EAS[!is.na(pheno_vad_EAS[,trait]),c(trait,"age","age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",paste0("p_value_",idx))]
    
    a <- as.formula(paste0(trait,"~p_value_",idx))
    roc_obj <- roc(a,data = d,quiet = TRUE,auc = TRUE)
    AUC <- auc(roc_obj)
    ci_result <- ci.auc(roc_obj,boot.n = 1000,method = "bootstrap",progress = "none")
    
    
    AUC.result <- data.frame(method = "CT_EAS",
                             AUC = AUC,
                             AUC_low = ci_result[1],
                             AUC_high = ci_result[3]
    )
    
    ## Save the AUC for the validation set w/ its confidence bounds, as well as the AUC tuning vector
    ct.result <- list(AUC.result,AUC_tun_vec)
  }
  ct.result_EAS <- AUC.result
  
  ct_result_rawAUC <- rbind(ct_result_rawAUC,rbind(ct.result_EUR,ct.result_NonEUR,ct.result_AFR,ct.result_EAS,ct.result_MIX,ct.result_SAS,ct.result_UNK))
  
}

rownames(ct_result_rawAUC) <- NULL
ct_result_rawAUC$Ancestry <- unlist(lapply(strsplit(ct_result_rawAUC$method,"_"),function(x) x[2]))
ct_result_rawAUC$Trait <- rep(c("Asthma","CAD","T2D","Breast","Prostate"),each = 7)
ct_result_rawAUC$method <- "CT_RawAUC"


a <- data.frame(Trait = NULL,Ancestry = NULL,Method = NULL,AUC = NULL,AUC_Low = NULL,AUC_High = NULL)

Trait <- "Asthma"
Ancestry <- "NonEur"

for(Trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  for(Ancestry in c("Eur","NonEur","UNK","MIX","AFR","SAS","EAS")){
    
    if(Ancestry == "Eur"){
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",Trait,"_CT_result_",toupper(Ancestry),".RData")) 
    }else{
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/CT/",Trait,"_CT_result_",Ancestry,".RData"))
    }
    if(sum(is.na(ct.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "CT",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "CT",AUC = ct.result[[1]]$AUC,AUC_Low = ct.result[[1]]$AUC_low,AUC_High = ct.result[[1]]$AUC_high))
    }
    
    if(Ancestry == "Eur"){
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",Trait,"_ldpred2_result_",toupper(Ancestry),".RData"))
    }else{
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LDPred2/",Trait,"_ldpred2_result_",Ancestry,".RData"))
    }
    if(sum(is.na(ldpred2.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LDPred2",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LDPred2",AUC = ldpred2.result[[1]]$AUC,AUC_Low = ldpred2.result[[1]]$AUC_low,AUC_High = ldpred2.result[[1]]$AUC_high))
    }
    
    if(Ancestry == "Eur"){
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",Trait,"_LASSOSUM2_result_",toupper(Ancestry),".RData"))
    }else{
      load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/LASSOSUM2/",Trait,"_LASSOSUM2_result_",Ancestry,".RData"))
    }
    if(sum(is.na(ldpred2.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LASSOSUM2",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "LASSOSUM2",AUC = LASSOSUM2.result[[1]]$AUC,AUC_Low = LASSOSUM2.result[[1]]$AUC_low,AUC_High = LASSOSUM2.result[[1]]$AUC_high))
    }
    
    
    load(paste0("/data/williamsjacr/UKB_WES_Phenotypes/Binary/Results/Combined_Common_PRS/",Trait,"_sl_result_All_",Ancestry,".RData"))
    if(sum(is.na(SL.result))>0){
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "SL_Common",AUC = NA,AUC_Low = NA,AUC_High = NA))
    }else{
      a <- rbind(a,data.frame(Trait = Trait,Ancestry = Ancestry,Method = "SL_Common",AUC = SL.result$AUC,AUC_Low = SL.result$AUC_low,AUC_High = SL.result$AUC_high))
    }
    
  }
}

a <- a[a$Method == "CT",]

ct_result_rawAUC <- as.data.frame(ct_result_rawAUC)
ct_result_rawAUC <- ct_result_rawAUC[,c("Trait","Ancestry","method","AUC","AUC_low","AUC_high")]
colnames(ct_result_rawAUC) <- c("Trait","Ancestry","Method","AUC","AUC_Low","AUC_High")
a <- rbind(a,ct_result_rawAUC)


results <- a
rm(a)

results$Ancestry <- toupper(results$Ancestry)
results$AUC <- as.numeric(results$AUC)

results <- results[!is.na(results$Ancestry),]

# a <- as.matrix(results)
# a <- as.data.frame(a)
# a$AUC <- as.numeric(a$AUC)
# a$AUC_Low <- as.numeric(a$AUC_Low)
# a$AUC_High <- as.numeric(a$AUC_High)

ggplot(results) +
  geom_bar(aes(x=Method, y=AUC,fill=Method), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Method, ymin=AUC_Low, ymax=AUC_High), width=0.4, colour="black", alpha=0.9) +  
  facet_grid(vars(Trait), vars(Ancestry)) + 
  ggtitle("Common Variants") + 
  #  ylim(0.45,1.05) + 
  ylab("AUC") + 
  theme_Publication() + 
  scale_fill_Publication()
