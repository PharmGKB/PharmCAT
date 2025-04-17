#!/bin/bash
#
#  Sets up PharmCAT in the current directory after checking for basic requirements.
#
#

set -e
set -u
set -o pipefail

PHARMCAT_VERSION=3.0.0
# Java 17 = 61 class major version
MIN_JAVA_VERSION=17
MIN_JAVA_CLASS_VERSION=61
MIN_PYTHON_VERSION=3.10
MIN_BGZIP_VERSION=1.18
MIN_BCFTOOLS_VERSION=1.18


if [[ $TERM && $TERM != 'dumb' ]]; then
  if command -v tput &>/dev/null; then
    GREEN=$(tput setaf 2; tput bold)
    YELLOW=$(tput setaf 3)
    RED=$(tput setaf 1)
    NORMAL=$(tput sgr0)
  fi
fi

function echo_red() {
  >&2 echo -e "$RED$*$NORMAL"
}

function echo_yellow() {
    >&2 echo -e "$YELLOW$*$NORMAL"
}

function echo_green() {
    echo -e "$GREEN$*$NORMAL"
}

darwin=false
case "$(uname)" in
    Darwin*)
        darwin=true
        ;;
esac


pipe_strip_newlines () {
  local str
  # read from stdin
  str="$(</dev/stdin)"
  str="${str//$'\r'/}"
  str="${str//$'\n'/}"
  echo "$str"
}

# compare function based on https://stackoverflow.com/a/4025065/1063501
# modified to use echo return values based on https://linuxize.com/post/bash-functions/#return-values
compare_versions () {
  if [[ $1 == "$2" ]]
  then
    echo "="
    return
  fi
  local IFS=.
  local i ver1=($1) ver2=($2)
  # fill empty fields in ver1 with zeros
  for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
  do
    ver1[i]=0
  done
  for ((i=0; i<${#ver1[@]}; i++))
  do
    if ((10#${ver1[i]:=0} > 10#${ver2[i]:=0}))
    then
      echo ">"
      return
    fi
    if ((10#${ver1[i]} < 10#${ver2[i]}))
    then
      echo "<"
      return
    fi
  done
  echo "="
  return
}


# check for pre-existing preprocessor directory
if [ -d "preprocessor" ]; then
    echo_yellow "Found pre-existing preprocessor subdirectory."
  	echo ""
    echo_yellow "Please delete this subdirectory and try again."
    echo ""
    exit 1
fi


echo "Looking for Java ${MIN_JAVA_VERSION} or greater..."
# xargs at the end removes newline
java_version=$(javap -verbose java.lang.String | grep "major version" | cut -d " " -f5 | pipe_strip_newlines)
if [[ $java_version -lt $MIN_JAVA_CLASS_VERSION ]]; then
	echo_red "Not found!"
	echo ""
	echo_yellow "Please install Java ${MIN_JAVA_VERSION} or greater on your system using your favorite package manager"
	echo_yellow "and try again."
	echo ""
	exit 1
fi

echo "Looking for Python 3.10 or greater..."
python_version=$(python3 -c "import platform; print(platform.python_version())")
op=$(compare_versions $python_version $MIN_PYTHON_VERSION)
if [[ $op = '<' ]]; then
	echo_red "Not found!"
	echo ""
	echo_yellow "Please install Python ${MIN_PYTHON_VERSION} or greater on your system using your favorite package manager"
	echo_yellow "and try again."
	echo ""
	exit 1
else
  echo "  Found Python ${python_version}"
fi

echo "Looking for bgzip ${MIN_BGZIP_VERSION} or greater..."
bgzip_version=$(bgzip --version | head -n 1 | awk '{print $3}')
op=$(compare_versions "${bgzip_version}" $MIN_BGZIP_VERSION)
if [[ $op = '<' ]]; then
	echo_red "Not found!"
	echo ""
	echo_yellow "Please install bgzip ${MIN_BGZIP_VERSION} or greater on your system and try again."
	echo ""
 	echo_yellow "bgzip can be found in htslib @ https://www.htslib.org/download/"
	if [[ $darwin == true ]]; then
	  echo_yellow "htslib is also available via homebrew"
	fi
	echo ""
	exit 1
else
  echo "  Found bgzip ${bgzip_version}"
fi

echo "Looking for bcftools ${MIN_BCFTOOLS_VERSION} or greater..."
bcftools_version=$(bcftools --version | head -n 1 | awk '{print $2}')
op=$(compare_versions "${bcftools_version}" $MIN_BCFTOOLS_VERSION)
if [[ $op = '<' ]]; then
	echo_red "Not found!"
	echo ""
	echo_yellow "Please install bcftools ${MIN_BCFTOOLS_VERSION} or greater on your system and try again."
	echo ""
	echo_yellow "bcftools can be found @ https://www.htslib.org/download/"
	if [[ $darwin == true ]]; then
	  echo_yellow "bcftools is also available via homebrew"
	fi
	echo ""
	exit 1
else
  echo "  Found bcftools ${bcftools_version}"
fi


RELEASE_URL=https://github.com/PharmGKB/PharmCAT/releases/download/v${PHARMCAT_VERSION}
PIPELINE_URL=${RELEASE_URL}/pharmcat-pipeline-${PHARMCAT_VERSION}.tar.gz

if command -v curl &>/dev/null; then
  GET="curl -fsSL '${PIPELINE_URL}' -o 'pipeline.tgz'"
elif command -v wget &>/dev/null; then
  GET="wget -nv '${PIPELINE_URL}' -O 'pipeline.tgz' >/dev/null 2>&1"
else
  echo_red "ERROR: Cannot find 'curl' nor 'wget' utility --  please install one of them"
	echo ""
  exit 1
fi


echo ""
echo "Downloading PharmCAT..."
set +e
eval $GET; status=$?
set -e
if [ $status -ne 0 ]; then
  echo ""
  echo_red "ERROR: Cannot download PharmCAT -- make sure you can connect to the internet"
  echo ""
  echo "Alternatively, you can try to download this file:"
  echo_yellow "  ${PIPELINE_URL}"
  echo "and un-tar it"
  echo ""
  exit 1
fi

tar -xzf pipeline.tgz
rm pipeline.tgz
chmod 755 pharmcat pharmcat_pipeline pharmcat_vcf_preprocessor

echo_green "Done!"
echo ""
echo "Don't forget to install PharmCAT VCF Preprocessor's Python dependencies."
echo ""
echo "You can do so with:"
echo_yellow "  pip3 install -r requirements.txt"
echo ""
echo "Or use your preferred Python virtual environment to do so."
echo ""
