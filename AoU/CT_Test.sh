%%writefile CT_Test.sh
#!/bin/bash

set -o errexit
set -o nounset

Rscript ${R_Script} ${BED_Full_File} ${BED_Ref_File} ${sumstats} ${all_phenotypes_file} ${all_train_file} ${all_tune_file} ${all_valid_file} ${anc} ${trait} ${OUTPUT_PATH}