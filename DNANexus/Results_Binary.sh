# load bigsnpr image

dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/CT/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/LDPred2_LASSOSum2/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/Combined_Common_PRS/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestRareVariantPRS/ -r
dx download UKB_PRS:/JW/UKB_Phenotypes/Results/Binary/BestPRS/ -r

mkdir Results_Binary

cp -a CT/. Results_Binary/
cp -a LDPred2_LASSOSum2/. Results_Binary/
cp -a Combined_Common_PRS/. Results_Binary/
cp -a BestRareVariantPRS/. Results_Binary/
cp -a BestPRS/. Results_Binary/

cd Results_Binary/
rm *.txt
rm *.log
rm *.sscore
cd ..

tar -czvf Results_Binary.tar.gz Results_Binary/

rm -r CT/
rm -r LDPred2_LASSOSum2/
rm -r Combined_Common_PRS/
rm -r BestRareVariantPRS/
rm -r BestPRS/