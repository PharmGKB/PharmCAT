#!/bin/bash
#
# Iterate through definition files and call test_gen.py to generate test VCF's.
# This script expects to be called within the PharmCAT repository (all paths are hard coded).
#
#
set -e
set -u
set -o pipefail


MISSING="false"

while getopts 'm' OPTION; do
  case "$OPTION" in
    m)
      MISSING="true"
      ;;

    ?)
      echo "script usage: $(basename $0) [-m]" >&2
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"



# cd to location of script
cd $(dirname $0)

if [ -z ${PHARMCAT_DATA_DIR+x} ]; then
  dataDir="../../../build"
else
  # expect PHARMCAT_DATA_DIR to be a relative directory
  dataDir="../../../${PHARMCAT_DATA_DIR}"
fi

DEFINITION_DIR="../../main/resources/org/pharmgkb/pharmcat/definition/alleles"
OUTPUT_DIR="${dataDir}/testVcf"

for file in "$DEFINITION_DIR"/*; do
  if [[ $file == *_translation.json ]]; then
    gene=$(echo $(basename $file) | cut -d "_" -f1)
    # TODO: remove after v1 is released
    # not performing testing on CYP2D6 for now
    if [[ $gene != "CYP2D6" ]]; then
      echo "$gene"
      mkdir -p ${OUTPUT_DIR}/${gene}
      if [[ $MISSING == "false" ]]; then
        ./test_gen.py $file ${OUTPUT_DIR}/${gene}
      else
        ./test_gen_missing.py $file ${OUTPUT_DIR}/${gene}
      fi
      if [ "$PHARMCAT_TEST_QUIET" != "true" ]; then
        echo ""
        echo ""
      fi
    fi
  fi
done
