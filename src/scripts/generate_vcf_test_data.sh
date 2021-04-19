#!/bin/bash
#
# Iterate through definition files and call test_gen.py to generate test VCF's.
# This script expects to be called within the PharmCAT repository (all paths are hard coded).
#
#

# cd to location of script
cd $(dirname $0)

DEFINITION_DIR="../main/resources/org/pharmgkb/pharmcat/definition/alleles"
VCF_DIR="../../testVcf"

for file in "$DEFINITION_DIR"/*; do
  if [[ $file == *_translation.json ]]; then
    gene=$(echo $(basename $file) | cut -d "_" -f1)
    echo "$gene"
    mkdir -p ${VCF_DIR}/${gene}
    python test_gen.py $file ${VCF_DIR}/${gene}
    echo ""
    echo ""
  fi
done
