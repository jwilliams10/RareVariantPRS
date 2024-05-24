%%writefile merge_bed.sh
#!/bin/bash

set -o errexit
set -o nounset

Rscript ${R_Script} ${INPUT_PATH} ${OUTPUT_PATH}