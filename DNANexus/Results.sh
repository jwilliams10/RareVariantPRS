# load bigsnpr image

dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/CT/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/LDPred2_LASSOSum2/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/Combined_Common_PRS/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestRareVariantPRS/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Continuous/BestPRS/ -r

mkdir Results

cp -a CT/. Results/
cp -a LDPred2_LASSOSum2/. Results/
cp -a Combined_Common_PRS/. Results/
cp -a BestRareVariantPRS/. Results/
cp -a BestPRS/. Results/

cd Results/
rm *.txt
rm *.sscore
cd ..

tar -czvf Results.tar.gz Results/

rm -r CT/
rm -r LDPred2_LASSOSum2/
rm -r Combined_Common_PRS/
rm -r BestRareVariantPRS/
rm -r BestPRS/