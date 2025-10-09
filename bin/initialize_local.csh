#!/usr/bin/tcsh

# install python packages
if ( ! -d "pcat_venv" ) then
  echo "Initializing venv (pcat_venv)"
  python3 -m venv pcat_venv
else
  echo "venv already initialized"
endif

source pcat_venv/bin/activate.csh

echo "Installing required python packages"
pip3 install -r preprocessor/requirements.txt
pip3 install -r preprocessor/tests/requirements.txt


# install samtools
set lsb_exists = `which lsb_release`
if ("$lsb_exists" != "") then
  set distro = `lsb_release -i -s`
  if ("$distro" == "Ubuntu") then
    echo "Installing required libraries..."
    sudo apt install libbz2-dev liblzma-dev libncurses5-dev libcurl4-gnutls-dev
  endif
endif
set called = ($_)
if ($#called > 1) then
  set script_path = $called[2]
else
  set script_path = $0
${script_path}/install_samtools.sh ~/bin
