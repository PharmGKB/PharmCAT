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
GENE=""

while getopts 'mg' OPTION; do
  case "$OPTION" in
    m)
      MISSING="true"
      ;;

    g)
      GENE=${2}
      ;;

    ?)
      echo "script usage: $(basename $0) [-m]" >&2
      exit 1
      ;;
  esac
done
shift "$((OPTIND -1))"



# cd to location of script
cd "$(dirname $0)"

if [ -z ${PHARMCAT_DATA_DIR+x} ]; then
  dataDir="../../../build"
else
  # expect PHARMCAT_DATA_DIR to be a relative directory (relative to PHARMCAT)
  dataDir="../../../${PHARMCAT_DATA_DIR}"
fi

DEFINITION_DIR="../../main/resources/org/pharmgkb/pharmcat/definition/alleles"
OUTPUT_DIR="${dataDir}/testVcf"
echo "Writing test data to ${OUTPUT_DIR}"

for file in "$DEFINITION_DIR"/*; do
  if [[ $file == *_translation.json ]]; then
    gene=$(basename "$file" | cut -d "_" -f1)
    if [[ -n $GENE ]]; then
      if [[ $GENE != "$gene" ]]; then
        continue
      fi
    fi
    echo "$gene"
    mkdir -p "${OUTPUT_DIR}/${gene}"
    if [[ $MISSING == "false" ]]; then
      ./test_gen.py "$file" "${OUTPUT_DIR}/${gene}"
    else
      ./test_gen_missing.py "$file" "${OUTPUT_DIR}/${gene}"
    fi
    if [ "${PHARMCAT_TEST_QUIET-defined}" != "true" ]; then
      echo ""
      echo ""
    fi
  fi
done
