%%writefile GWAS_SumStats.sh
#!/bin/bash

set -o errexit
set -o nounset

Rscript ${R_Script} ${BED_File} ${all_phenotypes_file} ${all_train_file} ${anc} ${trait} ${OUTPUT_PATH}