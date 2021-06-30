#!/bin/bash
#
# Iterate through definition files and call test_gen.py to generate test VCF's.
# This script expects to be called within the PharmCAT repository (all paths are hard coded).
#
#

# cd to location of script
cd $(dirname $0)

DEFINITION_DIR="../main/resources/org/pharmgkb/pharmcat/definition/alleles"
VCF_FILE="../../build/pharmcat_positions.vcf"
OUTPUT_DIR="../../build/testVcf"

for file in "$DEFINITION_DIR"/*; do
  if [[ $file == *_translation.json ]]; then
    gene=$(echo $(basename $file) | cut -d "_" -f1)
    # TODO: remove after v1 is released
    # not performing testing on CYP2D6 for now
    if [[ $gene != "CYP2D6" ]]; then
      echo "$gene"
      mkdir -p ${OUTPUT_DIR}/${gene}
      ./test_gen.py $file ${VCF_FILE} ${OUTPUT_DIR}/${gene}
      echo ""
      echo ""
    fi
  fi
done
