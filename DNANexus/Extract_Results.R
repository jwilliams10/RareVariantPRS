# load bigsnpr image

system("mkdir Results")

setwd("Results/")

for(trait in c("BMI","HDL","LDL","Height","TC","logTG")){
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/CT/",trait,"Best_Betas.csv"))
  system(paste0("mv ",trait,"Best_Betas.csv ",trait,"_Best_Betas_CT.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/",trait,"Best_Betas.csv"))
  system(paste0("mv ",trait,"Best_Betas.csv ",trait,"_Best_Betas_CV_SL.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/",trait,"_Coding_Best_Betas.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/",trait,"_Noncoding_Best_Betas.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/",trait,"_Best_Betas.csv"))
  system(paste0("mv ",trait,"_Best_Betas.csv ",trait,"_Best_Betas_RV_SL.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestPRS/",trait,"Best_Betas.csv"))
  system(paste0("mv ",trait,"Best_Betas.csv ",trait,"_Best_Betas_CV_RV.csv"))
}

setwd("~")
system("tar -czvf Results.tar.gz Results/")

system("dx upload Results.tar.gz --path UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/")







# load bigsnpr image

system("mkdir Results_Binary")

setwd("Results_Binary/")

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/CT/",trait,"Best_Betas.csv"))
  system(paste0("mv ",trait,"Best_Betas.csv ",trait,"_Best_Betas_CT.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LASSOSum.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/",trait,"Best_Betas_LDPred2.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/",trait,"Best_Betas.csv"))
  system(paste0("mv ",trait,"Best_Betas.csv ",trait,"_Best_Betas_CV_SL.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/",trait,"_Coding_Best_Betas.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/",trait,"_Noncoding_Best_Betas.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/",trait,"_Best_Betas.csv"))
  system(paste0("mv ",trait,"_Best_Betas.csv ",trait,"_Best_Betas_RV_SL.csv"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestPRS/",trait,"Best_Betas.csv"))
  system(paste0("mv ",trait,"Best_Betas.csv ",trait,"_Best_Betas_CV_RV.csv"))
}

setwd("~")
system("tar -czvf Results_Binary.tar.gz Results_Binary/")

system("dx upload Results_Binary.tar.gz --path UKB_PRS:JW/UKB_Phenotypes/Results/Binary/")







system("mkdir Binary_PRS_Validation")

setwd("Binary_PRS_Validation/")

for(trait in c("Asthma","CAD","T2D","Breast","Prostate")){
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/",trait,"_BestPRS.csv"))
}

setwd("~")
system("tar -czvf Binary_PRS_Validation.tar.gz Binary_PRS_Validation/")

system("dx upload Binary_PRS_Validation.tar.gz --path UKB_PRS:JW/UKB_Phenotypes/Results/Binary/")








system("mkdir Continuous_PRS_Validation")

setwd("Continuous_PRS_Validation/")

for(trait in c("BMI","HDL","LDL","Height","TC","logTG")){
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/",trait,"_Best_Validation_All.txt"))
  system(paste0("dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/",trait,"_BestPRS.csv"))
}

setwd("~")
system("tar -czvf Continuous_PRS_Validation.tar.gz Continuous_PRS_Validation/")

system("dx upload Continuous_PRS_Validation.tar.gz --path UKB_PRS:JW/UKB_Phenotypes/Results/Continuous/")


