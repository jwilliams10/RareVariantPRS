#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=144:00:00
#SBATCH --array=22
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.0

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/gds_processing/vcf_to_gds.R ${SLURM_ARRAY_TASK_ID} > vcf_to_gds"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/gds_processing/Add_QC_label.R > Add_QC_label.Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/gds_processing/Varinfo_gds.R ${SLURM_ARRAY_TASK_ID} > Varinfo_gds"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/gds_processing/Annotate.R ${SLURM_ARRAY_TASK_ID} > Annotate"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/gds_processing/gds2agds.R ${SLURM_ARRAY_TASK_ID} > gds2agds"${SLURM_ARRAY_TASK_ID}".Rout
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/RareVariantPRS/Simulation_Study/FullData/Extract_chr22/gds_processing/Association_Analysis_Prestep.R > Association_Analysis_Prestep.Rout