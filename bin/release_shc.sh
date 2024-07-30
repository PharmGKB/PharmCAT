 #!/bin/bash
 #
 # Performs a "release" for SHC PharmCAT.
 #
 #
 set -e
 set -u
 set -o pipefail


pcat_version=$(git tag -l --sort=-creatordate | head -n 1)
shc_version=${pcat_version}-shc

echo "Releasing ${shc_version}"
git checkout -b shc-${pcat_version}

echo ""
echo "Updating data"
make updateData
git commit -m "chore: add base data"

echo ""
echo "Subsetting data"
make updateShc
git add src/main/resources/ pharmcat_positions.*

echo ""
echo "Updating version in preprocessor"
sed -i -E "s/PHARMCAT_VERSION = '.+'/PHARMCAT_VERSION = '${shc_version}'/" preprocessor/preprocessor/common.py
git add preprocessor/preprocessor/common.py

git commit -m "feat(shc): subset data"


echo ""
echo "Tagging ${shc_version}"
git tag ${shc_version}

echo ""
echo "Done."
echo ""


