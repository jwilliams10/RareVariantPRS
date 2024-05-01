# load bigsnpr image

mkdir Results

cd Results/

for trait in BMI HDL LDL Height TC logTG;
do
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/CT/${trait}Best_Betas.csv
  mv ${trait}Best_Betas.csv ${trait}_Best_Betas_CT.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}Best_Betas_LASSOSum.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/${trait}Best_Betas_LDPred2.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/${trait}Best_Betas.csv
  mv ${trait}Best_Betas.csv ${trait}_Best_Betas_CV_SL.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}_Coding_Best_Betas.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}_Noncoding_Best_Betas.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/${trait}_Best_Betas.csv
  mv ${trait}_Best_Betas.csv ${trait}_Best_Betas_RV_SL.csv
  dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestPRS/${trait}Best_Betas.csv
  mv ${trait}Best_Betas.csv ${trait}_Best_Betas_CV_RV.csv
done

cd ..

tar -czvf Results.tar.gz Results/
