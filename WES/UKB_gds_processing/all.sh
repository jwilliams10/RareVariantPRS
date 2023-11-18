#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=1-22
#SBATCH --mem-per-cpu=40G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/UKB_gds_processing/vcf_to_gds.R ${SLURM_ARRAY_TASK_ID} > vcf_to_gds"${SLURM_ARRAY_TASK_ID}".Rout 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/UKB_gds_processing/Add_QC_label.R ${SLURM_ARRAY_TASK_ID} > Add_QC_label"${SLURM_ARRAY_TASK_ID}".Rout 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/UKB_gds_processing/Varinfo_gds.R ${SLURM_ARRAY_TASK_ID} > Varinfo_gds"${SLURM_ARRAY_TASK_ID}".Rout 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/UKB_gds_processing/Annotate.R ${SLURM_ARRAY_TASK_ID} > Annotate"${SLURM_ARRAY_TASK_ID}".Rout 
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/WES/UKB_gds_processing/gds2agds.R ${SLURM_ARRAY_TASK_ID} > gds2agds"${SLURM_ARRAY_TASK_ID}".Rout 
