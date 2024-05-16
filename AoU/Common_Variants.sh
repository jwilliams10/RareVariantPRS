%%writefile Common_Variants.sh
#!/bin/bash

set -o errexit
set -o nounset

Rscript ${R_Script} ${CHR} ${BED_File} ${BIM_File} ${FAM_File} ${Ancestry_File} ${OUTPUT_PATH}