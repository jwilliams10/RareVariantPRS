for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; 
do
  dx run app-swiss-army-knife -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.fam -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bed -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bim -icmd="plink -bfile ukb24305_c${chrom}_b0_v1 --maf 0.01 --hwe 0.000001 --geno 0.02 --mind 0.05 --make-bed --out chr${chrom}_filtered_common" -y --destination UKB_PRS:JW/Clean_Data/ --instance-type mem1_ssd1_v2_x72
done

for chrom in 1 2 3 4 5; 
do
  dx run app-swiss-army-knife -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.fam -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bed -iin=UKB_PRS:Bulk/'Whole genome sequences'/'Population level WGS variants, PLINK format - interim 200k release'/ukb24305_c${chrom}_b0_v1.bim -icmd="plink -bfile ukb24305_c${chrom}_b0_v1 --maf 0.01 --hwe 0.000001 --geno 0.02 --mind 0.05 --make-bed --out chr${chrom}_filtered_common" -y --destination UKB_PRS:JW/Clean_Data/ --instance-type mem1_hdd1_v2_x36
done